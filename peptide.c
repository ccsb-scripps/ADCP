/*
** This contains the description of polypeptide structure
** and the routines for IO in PDB-like format.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<float.h>

#include"error.h"
#include"params.h"
#include"aadict.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"
#include"energy.h"

#define FAR 14887843.138588

/***********************************************************/
/****                     CONSTANTS                     ****/
/***********************************************************/

/*
  Polypeptide structure (Engh and Huber, 1991, 2001). !!! DO NOT EDIT !!!
*/

/* fixed Ca-Cb bond length and Ca-Ca distance */
const double cacb = 1.532, caca = 3.819, caca_ = 2.922322706341652474;

/* fixed trans peptide bond atom coordinates in local frames */
vector ca = { 0.0, 0.0, 3.819 };
vector h = { 0.0, 1.357, -1.643 };
vector n = { 0.0, 0.388, -1.406 };
vector c = { 0.0, -0.520, 1.434 };
vector o = { 0.0, -1.730, 1.650 };

/* fixed cis (proline) peptide bond atom coordinates */
vector ca_ = { 0.0, 1.849, 2.263 };
vector h_ = { 0.0, -1.812, 1.083 };
vector n_ = { 0.0, -1.460, 0.150 };

/* normalized bond vectors */
vector n1, c1, n1_, cn1, cn1_;

/* approximate reciprocal Ca-N and Ca-C bond lengths and skew between them */
double ican = 0.6854;
double icac = 0.6556;
double ican_ = 0.6813;
double skew = -3.063;
double skew_ = -1.121;


/***********************************************************/
/****            INITIALISATION OF CONSTANTS            ****/
/***********************************************************/


/* assure the correspondence between the coordinates, lengths, and angles */
//void peptide_init(char *prm)
void peptide_init()
{
	vector a;

	ican = 1.0 / sqrt(square(n));
	icac = 1.0 / sqrt(square(c));
	ican_ = 1.0 / sqrt(square(n_));

	scale(n1, ican, n);
	scale(c1, icac, c);
	scale(n1_, ican_, n_);

	skew = atan2(n[1] * c[2] - n[2] * c[1], n[1] * c[1] + n[2] * c[2]);
	skew_ = atan2(n_[1] * c[2] - n_[2] * c[1], n_[1] * c[1] + n_[2] * c[2]);

	/* printf("%g %g %g\n", ican, icac, skew); */

	add(a, ca, n);
	subtract(a, a, c);
	normalize(a);
	castvec(cn1, a);

	add(a, ca_, n_);
	subtract(a, a, c);
	normalize(a);
	castvec(cn1_, a);
}


/***********************************************************/
/****              PEPTIDE GEOMETRY  TESTS              ****/
/***********************************************************/


/* Detect disulfide bridges. */
void chkssbond( AA *a, int count)
{
	int i, j;
	double p;

	for (i = 1; i < count; i++) {
		if (a[i].id != 'C')
			continue;

		for (j = 1; j < i; j++)
			if (a[j].id == 'C' &&
			    (p = distance(a[j].g, a[i].g)) < 6.0)
				fprintf(stderr, "Disulfide bond C%d-C%d : %f\n",
					a[j].num, a[i].num, sqrt(p));
	}
}

/* Detect missing atoms in an amino acid, by comparing the AA.etc
   with the list of atoms that should be present in the amino acid.
   CAUTION! This will only work, if the AA.etc is up to date. */
static void check_missing_atoms(const AA a, model_params *mod_params) {

	int mask = N__ | CA_ | C__ | O__ ;

	/* beta carbon */
	if (a.id != 'G')
		mask |= CB_ ;
	if (a.id != 'P')
		mask |= H__ ;

	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		/* gamma atom(s) */
#ifdef LINUS_1995
		/* LINUS 1995 doesn't have CG for PRO */
		if (a.id != 'G' && a.id != 'A' && a.id != 'P')
#else
		if (a.id != 'G' && a.id != 'A')
#endif
			mask |= G__;
		if (a.id == 'I' || a.id == 'T' || a.id == 'V')
			mask |= G2_;

	}
	/* CAUTION! Because gamma atoms have been added, this is not backwards compatible!
	   When reading an old file, this will always complain about missing gamma atoms. */
	if (mask & ~a.etc) {
		fprintf(stderr, "Partial residue  %s %5d : %#x\n", aa123(a.id), a.num, mask & ~a.etc);
	}
}

/* Calculate the chirality of an amino acid around the CA atom.
   If the d(CA-CB) > 2 Angstroms, report that this is an odd residue.
   Levo: (NCA.CCAxCBCA) ~ 2.5 Angstroms, but accept -1.0 < (NCA.CCAxCBCA) < 4.0.
   Dextro: accept -4.0 < (NCA.CCAxCBCA) < -1.0.
   Set the chirality in AA.etc, and report if it is not a levo residue. */
static void calculate_aa_chirality(AA *a) {

	a->etc |= LEV;	/* initialization */

	if (a->id == 'G')
		return;

	vector x, y, z;
	subtract(x, a->n,  a->ca);
	subtract(y, a->c,  a->ca);
	subtract(z, a->cb, a->ca);
	double p = square(z);	/* if less than 4, Cb is reliable */
	double q = triprod(x, y, z);	/* about 2.5 for levo amino acids */

	if (p > 4.0 || q < -4.0 || 4.0 < q)
		fprintf(stderr, "Odd amino acid   %s %5d : %g %g\n",
			aa123(a->id), a->num, sqrt(p), q);

	/* dextro amino acids have negative q of about -2.5 */
	if (-4.0 < q && q < -1.0 && p < 4.0) {
		fprintf(stderr, "Dextro residue   %s %5d : %g %g\n",
			aa123(a->id), a->num, sqrt(p), q);
		a->etc &= ~LEV;
	}
}

/* Check for gaps in the chain, that is if d(CA_i,CA_i+1) > 20 Angstroms
   or d(C_i,N_i+1) > 4 Angstroms.
   Also calculate whether the peptide bond connecting an amino acid
   to the previous amino acid is cis/trans and save it in AA.etc. 
   Ignore the false connections through chain breaks marked by chainid. */
static void check_chain_connectivity(AA *a, int count) {

	double p, q;
	a[1].etc &= ~CIS;
	for (int i = 2; i < count; i++) {

		/* ignore amino acid pairs in different chains for multi-chain proteins */
		if (a[i].chainid != a[i - 1].chainid) {
			fprintf(stderr,"Skipping chain break between amino acids %d and %d.\n", i-1, i);
			/* set default trans peptide bond and skip */
			a[i].etc &= ~CIS;
			continue;
		}

		p = distance(a[i].n, a[i - 1].c);
		q = distance(a[i].ca, a[i - 1].ca);

		if (p > 4.0 || q > 20.0) {
			fprintf(stderr, "Chain gap %4d   %s %5d : %g %g\n",
				a[i].num - a[i - 1].num - 1, aa123(a[i].id),
				a[i].num, sqrt(p), sqrt(q));
			fprintf(stderr, "      %d %d %s // %d %d %s\n", i-1, a[i-1].num, aa123(a[i-1].id), i, a[i].num, aa123(a[i].id));
		}
		/* cis-peptides have smaller Ca-Ca of 2.922 A */
		if (q < 11.0 && a[i].id != 'P')
			fprintf(stderr, "Cis peptide bond %s %5d : %g %g\n",
				aa123(a[i].id), a[i].num, sqrt(p), sqrt(q));

		if (q < 11.0)
			a[i].etc |= CIS;
		else
			a[i].etc &= ~CIS;
	}

}

/* Invalidate nonexisting atoms of an amino acid by setting the coordinates
   of those atoms to FAR. */
static void invalidate_nonexisting_atoms(AA *a, model_params *mod_params) {

	/* GLY: no beta carbon */
	if (a->id == 'G')
		scale(a->cb, FAR * a->num, h);

	/* PRO: no backbone hydrogen */
	if (a->id == 'P')
		scale(a->h, -FAR * a->num, h);

	/* gamma atom(s) */
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
#ifdef LINUS_1995
		/* LINUS 1995 doesn't have CG for PRO */
		if (a->id == 'G' || a->id == 'A' || a->id == 'P')
#else
		if (a->id == 'G' || a->id == 'A')
#endif
			scale(a->g, FAR * a->num, ca);
		if (a->id != 'I' && a->id != 'T' && a->id != 'V')
			scale(a->g2, FAR * a->num, ca);
	} else {
		scale(a->g, FAR * a->num, ca);
		scale(a->g2, FAR * a->num, ca);
	}
}

/* Mark hydrophobic and amphipathic amino acids in AA.etc, according to
   R. Srinivasan et al., PNAS 96(25), 14258--14263 (1999).
   CAUTION! The hydrophobic atoms within amino acids are set
   by initialize_one_sidechain_properties in aadict.c. */
static void set_aa_hydrophobicity(AA *a) {

	a->etc &= ~HYDROPHOBIC;
	a->etc &= ~AMPHIPATHIC;
	if ( /* hydrophobic amino acids */
		a->id == 'C' || a->id == 'I' ||
		a->id == 'L' || a->id == 'M' ||
		a->id == 'F' || a->id == 'W' ||
		a->id == 'V')// || a->id == 'P')
		a->etc |= HYDROPHOBIC;
	if ( /* amphipathic amino acids */
		a->id == 'A' || a->id == 'H' ||
		a->id == 'T' || a->id == 'Y' || a->id == 'P' )
		a->etc |= AMPHIPATHIC;
}

/* Mark electrostatic amino acids in AA.etc, according to
   the total charge of their side chains.
   CAUTION! The electrostatic atoms within amino acids are set
   by initialize_one_sidechain_properties in aadict.c. */
static void set_aa_electrostaticity(AA *a) {

	a->etc &= ~ELECTROSTATIC;
	if ( /* charged amino acids */
		a->id == 'D' || a->id == 'E' ||
		a->id == 'K' || a->id == 'R' ) {
		a->etc |= ELECTROSTATIC;
	}
}

/* Detect irregularities in polypeptide chain and set various bits of information
   in AA.etc.  Only the invalid atoms' coordinates are touched. */
void chkpeptide(AA *a, int count, model_params *mod_params)
{
	int i;
	double p;

	/* detect missing atoms */
	for (i = 1; i < count; i++) check_missing_atoms(a[i], mod_params);

	/* chain chirality check */
	for (i = 1; i < count; i++) calculate_aa_chirality(&(a[i]));

	/* chain connectivity check */
	check_chain_connectivity(a, count);

	/* psi dihedral check */
	for (i = 1; i < count - 1; i++) {
		/* skip multi-chain proteins */
		p = dihedral_4(a[i].n, a[i].ca, a[i].c, a[i + 1].n);
		if (-2.0 * M_PI_3 < p && p < M_PI_3)
			a[i].etc &= ~PSI;
		else
			a[i].etc |= PSI;
	}
	a[count-1].etc &= ~PSI;

	/* disulfide bonds */
	chkssbond(a, count);

	/* re-enumerate and invalidate meaningless atoms */
	for (i = 1; i < count; i++) {
		a[i].num = i;
		invalidate_nonexisting_atoms(&(a[i]), mod_params);
	}

	/* mark hydrophobic and amphipathic, and electrostatic amino acids */
	/* hydrophobic and electrostatic atoms within amino acids are set by initialize_one_sidechain_properties in aadict.c */
	for (i = 1; i < count; i++) {
		set_aa_hydrophobicity(&(a[i]));
		set_aa_electrostaticity(&(a[i]));
	}

}


/***********************************************************/
/****          AMINO ACID AND PEPTIDE BUILDING          ****/
/****             AND PEPTIDE  MODIFICATION             ****/
/***********************************************************/

/* Calculate the CA position of amino acid *a*, given amino acid *b* and
   the CA_a->CA_b vector.  Amino acid *a* follows *b* in the sequence. */
void carbonate_f(AA *a, AA *b, triplet x)
{
	if (a->etc & CIS) {
		fling(a->ca, b->ca, ca_[1], x[1]);
		fling(a->ca, a->ca, ca_[2], x[2]);
	} else
		fling(a->ca, b->ca, caca, x[2]);
}

/* Calculate the CA position of amino acid *a*, given amino acid *b* and
   the CA_a->CA_b vector.  Amino acid *a* precedes *b* in the sequence. */
void carbonate_b(AA *a, AA *b, triplet x)
{
	if (b->etc & CIS) {
		fling(a->ca, b->ca, -ca_[1], x[1]);
		fling(a->ca, a->ca, -ca_[2], x[2]);
	} else
		fling(a->ca, b->ca, -caca, x[2]);
}

/* Calculate the position of the CB atom, given the CA->N and CA->C vectors.
   The position of the CB atom depends on the tau = <(N-CA-C) angle.
   TODO: needs checking these numbers. */
static double *branch(vector zz, vector can, vector cac, int levo)
{
	double p, q;
	vector xy;

	/* normalize bond vectors */
	scale(can, ican, can);
	scale(cac, icac, cac);

    add(zz, can, cac);
	crossprod(xy, can, cac);
	q = dotprod(can, cac);

	/* -CA- angles are flexible synchronous scissors */
	p = -.54835422911833237471;
	q = .56844147803928796138 - .75614562856111715263 * q;

	/* levo-dextro isomerism alters one sign */
	if (!levo)
		q = -q;

	lincomb(zz, p * cacb, zz, q * cacb, xy);

	return zz;
}

/* Calculate and set the gamma atom coordinates to the amino acid, given *chi*, that is
   the N-CA-C-G dihedral angle, and the side chain properties library.
   CAUTION!  The AA.etc mask whether the CB atom exists will not be overwritten here. */
void gammalate(AA *a, double chi, int which_gamma, model_params *mod_params)
{
	double c1, r, theta;
	triplet x;
	vector b;
	char routine_name[]="gammalate";
	char error_string[DEFAULT_LONG_STRING_LENGTH]="";

	if (chi == DBL_MAX) return;

	/* set the CB-G distance and CA-CB-G angle from the dictionary */
	if (beta_gamma_dist(a->id,which_gamma,&r,&theta,mod_params->sidechain_properties) != 0) {
	   sprintf(error_string,"%s: could not calculate r and theta.",routine_name);
	   stop(error_string);
	}

	/* local frame */
	subtract(x[2], a->cb, a->ca);
	//fprintf(stderr,"cb - ca = [%g %g %g] - [%g %g %g]\n",a->cb[0],a->cb[1],a->cb[2],a->ca[0],a->ca[1],a->ca[2]);
	normalize(x[2]);
	subtract(x[0], a->n, a->ca);
	c1 = dotprod(x[0], x[2]);
	fling(x[0], x[0], -c1, x[2]);
	normalize(x[0]);
	crossprod(x[1], x[2], x[0]);

	/* spherical coordinates */
	/* in the Euler coordintes PI-(ca-cb-g)angle is the theta angle */
	sphereframe(b, x, r, M_PI-theta, chi);
	if (which_gamma==1) {
	   add(a->g, a->cb, b);
	} else if (which_gamma==2) {
	   add(a->g2, a->cb, b);
	} else {
	   sprintf(error_string,"%s: invalid value for which gamma (%d).",routine_name,which_gamma);
	   stop(error_string);
	}
}


/* set backbone atoms and side chain (Cb) */
void acidate(AA *a, triplet xnt, triplet xct, simulation_params *sim_params)
{
	double *nn, *hh;
	vector zz, can, cac;	/* bonds Ca-N and Ca-C */
	//double chi1, chi2;

	if (a->etc & CIS) {
		hh = h_;
		nn = n_;
	} else {
		hh = h;
		nn = n;
	}

	/* backbone atoms */
	if (a->id != 'P') {
		vectortriplet(zz, hh, xnt);
		add(a->h, a->ca, zz);
	} else
		scale(a->h, -FAR * a->num, h);	/* meaninglessly far */

	vectortriplet(can, nn, xnt);
	add(a->n, a->ca, can);

	vectortriplet(zz, o, xct);
	add(a->o, a->ca, zz);

	vectortriplet(cac, c, xct);
	add(a->c, a->ca, cac);

	/* side chain */
	/* initialise meaninglessly far all */
	scale(a->cb, FAR * a->num, h);
	if ((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) {
		scale(a->g, FAR * a->num, h);
		scale(a->g2, FAR * a->num, h);
	}
	if (a->id != 'G') {
		/* beta carbon */
		branch(zz, can, cac, a->etc & LEV);
		add(a->cb, a->ca, zz);
		//fprintf(stderr,"zz = [%g %g %g]\n",zz[0],zz[1],zz[2]);
		//fprintf(stderr,"a->ca = [%g %g %g]\n",a->ca[0],a->ca[1],a->ca[2]);
		//fprintf(stderr,"a->cb = [%g %g %g]\n",a->cb[0],a->cb[1],a->cb[2]);
		/* possibly gamma atoms */
		if ((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) {
#ifdef LINUS_1995
			/* LINUS 1995 doesn't have CG for PRO */
			if (a->id != 'A' && a->id != 'P' ) {
#else
			if (a->id != 'A') {
#endif
				/* gamma atom */
				//fprintf(stderr,"gammalate, %d %c\n",a->num,a->id);
				//fprintf(stderr,"a->cb = [%g %g %g]\n",a->cb[0],a->cb[1],a->cb[2]);
				gammalate(a, a->chi1, 1, &(sim_params->protein_model));
				if (a->id == 'I' || a->id == 'T' || a->id == 'V' ) {
					/* second gamma atom */
					gammalate(a,a->chi2,2, &(sim_params->protein_model));
				}
			}
		}
	}

	/* fprintf(stderr, "%g %g\n",
	   cosine(zz, xaa[pos - 1][2]), cosine(zz, xaa[pos][2])); */
}

static double *atomate(vector a, vector q, vector b, vector c, vector o, model_params *mod_params)
{
	lincomb(a, q[0], b, q[1], c);
	fling(a, a, q[2], o);

	return a;
}

/* Fill in missing backbone H, O, N and C atoms, the existing CA, H, O, N, C atoms are not touched.
   Map the CB atoms onto the CRANKITE model.
   Report if there are any missing gamma atoms, but do not touch them (yet). */
/* TODO: also map gamma atoms (warning for now). */
void fixpeptide(AA *a, int count, model_params *mod_params)
{
	int i;
	double q;
	vector x, y, qq;
	
    //fprintf(stderr,"FIXING PEPTIDE\n");

	/* occasional nitrogen, needs CA and C */
	for (i = 2; i < count; i++) {
		/* skip chain breaks of multi-chain proteins */
		if (a[i].chainid != a[i-1].chainid)
			continue;

		if (a[i].etc & N__ || ~a[i - 1].etc & C__)
			continue;

		q = distance(a[i].ca, a[i - 1].ca);
		if (q < 11.0) { /* cis peptide bond */
			subtract(x, c, ca_);
			scale(y, -1.0, ca_);
			twocomps(qq, n_, y, x);
		} else if (q < 20.0) { /*trans peptide bond */
			subtract(x, c, ca);
			scale(y, -1.0, ca);
			twocomps(qq, n, y, x);
		} else
			continue;

		atomate(a[i].n, qq, a[i - 1].ca, a[i - 1].c, a[i].ca, mod_params);
		a[i].etc |= N__;
	}

	/* occasional carbon, needs CA and N */
	for (i = 1; i < count - 1; i++) {
		/* skip chain breaks of multi-chain proteins */
		if (a[i].chainid != a[i+1].chainid)
			continue;

		if (a[i].etc & C__ || ~a[i + 1].etc & N__)
			continue;

		q = distance(a[i].ca, a[i + 1].ca);
		if (q < 11.0) {
			add(x, n_, ca_);
			twocomps(qq, c, ca_, x);
		} else if (q < 20.0) {
			add(x, n, ca);
			twocomps(qq, c, ca, x);
		} else
			continue;

		atomate(a[i].c, qq, a[i + 1].ca, a[i + 1].n, a[i].ca, mod_params);
		a[i].etc |= C__;
	}

	/* hydrogen */
	/* works with multi-chain proteins */
	for (i = 1; i < count; i++) {
		/* starting hydrogens of all chains */
		/* WARNING: Could have more than 1 positions, we will just add it along the elongated backbone */
		if (a[i].chainid != a[i-1].chainid || i==1) {
			if (~(a[i].etc & H__) && a[i].id != 'P') {
				subtract(x, a[i].ca, a[i].c);
				fling(a[i].h,a[i].n, sqrt(distance(h,n)/square(x)), x);
				a[i].etc |= H__;
			}
		/* other hydrogens */
		} else {
			subtract(x, c, ca);
			twocomps(qq, h, n, x);
			if (a[i].etc & H__ || a[i].id == 'P')
				continue;
			atomate(a[i].h, qq, a[i].n, a[i - 1].c, a[i].ca, mod_params);
			a[i].etc |= H__;
		}
	}

	/* oxygen */
	add(x, n, ca);
	twocomps(qq, o, c, x);
	for (i = 1; i < count - 1; i++) {
		if (a[i].etc & O__)
			continue;
		if (a[i].chainid != a[i+1].chainid) continue;

		atomate(a[i].o, qq, a[i].c, a[i + 1].n, a[i].ca, mod_params);
		a[i].etc |= O__;
	}

	/* beta-carbon */
	for (i = 1; i < count; i++) {
		if (a[i].etc & CB_ || a[i].id == 'G' || ~a[i].etc & (N__ | C__))
			continue;
		if (a[i].etc & G__ || a[i].etc & G2_)
			fprintf(stderr,"WARNING: existing gamma atom(s), but missing beta carbon of residue %c%d\n",a[i].id,i);
		subtract(x, a[i].n, a[i].ca);
		subtract(y, a[i].c, a[i].ca);
		branch(qq, x, y, LEV);
		add(a[i].cb, a[i].ca, qq);
		a[i].etc |= CB_;
	}

	/* gamma atoms */
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		/* TODO: gamma atoms (warning for now) */
		for (i = 1; i < count; i++) {
#ifdef LINUS_1995
			/* LINUS 1995 doesn't have CG for PRO */
			if (a[i].etc & G__ || a[i].id == 'G' || a[i].id == 'A' || a[i].id == 'P' || ~a[i].etc & (N__ | C__))
#else
			if (a[i].etc & G__ || a[i].id == 'G' || a[i].id == 'A' || ~a[i].etc & (N__ | C__))
#endif
				continue;
			fprintf(stderr,"WARNING: missing gamma atom of residue %c%d\n",a[i].id,i);
		}
		for (i = 1; i < count; i++) {
			if (a[i].etc & G2_ || ( a[i].id != 'I' && a[i].id != 'V' && a[i].id != 'T' ) || ~a[i].etc & (N__ | C__))
				continue;
			fprintf(stderr,"WARNING: missing gamma atom of residue %c%d\n",a[i].id,i);
		}
	}

}

/* Adjust chi angle to the closest angle among 180, 60 and -60 degrees,
   the 3 stable sidechain dihedral angles */
void adjust_to_closest_state(double *chi, char aa_type) {
if(aa_type != 'P'){
	while (*chi < (-2.0*M_PI_3)) *chi +=2.0*M_PI;
	while (*chi > (4.0*M_PI_3)) *chi -=2.0*M_PI;

	if (*chi < (-2.0*M_PI_3)) {
	    stop("adjust_to_closest_state 1");
	} else if (*chi <= 0.0) {
	    *chi = -M_PI_3;
	} else if (*chi <= 2.0*M_PI_3) {
	    *chi = M_PI_3;
	} else if (*chi <= 4.0*M_PI_3) {
	    *chi = M_PI;
	} else {
	    stop("adjust_to_closest_state 2");
	}
      }
    else{
		while (*chi < (-1.0 * M_PI) ) *chi +=2.0*M_PI;
		while (*chi > M_PI ) *chi -=2.0*M_PI;
	    if( *chi < 0 ) *chi = -M_PI_6;
		else *chi = M_PI_6;
		
	}
}

/* Initialising the sidechain dihedral angle(s) of an amino acid,
   and saving them in the amino acid structure.  This will enable
   us to regenerate the side chain with the right dihedral angle
   at any time when the backbone moves, without having to keep saving them.
   The possibilities are:
     calculate the original dihedral angle(s), and save it as it is, or
     calculate the original dihedral angle(s), and save the closest value of +/-60 or 180 degrees, or
     use the default, most common value for the amino acid.
   Which one is used is set in mod_params. */
void initialise_sidechain_dihedral_angles(AA *a, model_params *mod_params ) {

	double chi1 = DBL_MAX, chi2 = DBL_MAX; /* set to DBL_MAX */
	char error_string[DEFAULT_LONG_STRING_LENGTH]="";
#ifdef LINUS_1995
	/* LINUS 1995 doesn't have CG for PRO */
	if (a->id != 'G' && a->id != 'A' && a->id != 'P' ) {
#else
	if (a->id != 'G' && a->id != 'A') {
#endif
		/* gamma atom */
		if (mod_params->use_original_gamma_atoms && (a->etc & G__) && (a->etc & CB_) && (a->etc & CA_) && (a->etc & N__)) {
			/* calc original chi, if possible */
			chi1 = dihedral_4(a->n,a->ca,a->cb,a->g);
		} 
		else { /* otherwise, use default most common value */
	//chi1 = sidechain_dihedral(a->id,1, mod_params->sidechain_properties);
    /*Pick random starting side chain dihedral */
      chi1 = sidechain_dihedral(a->id, mod_params->sidechain_properties);
		}
		if (mod_params->use_3_states) adjust_to_closest_state(&chi1,a->id);
		a->chi1 = chi1;
		/* second gamma atom */
		if (a->id == 'I' || a->id == 'T' || a->id == 'V' ) {
			if (mod_params->use_original_gamma_atoms && (a->etc & G2_) && (a->etc & CB_) && (a->etc & CA_) && (a->etc & N__)) {
				/* calc original chi2, if possible */
				chi2 = dihedral_4(a->n,a->ca,a->cb,a->g2);
			} 
			else { /* otherwise, use default most common value */
				//chi2 = sidechain_dihedral(a->id,2, mod_params->sidechain_properties);
				/* set chi2 based on chi1 */
				chi2 = sidechain_dihedral2(a->id,chi1, mod_params->sidechain_properties);
			}
			if (mod_params->use_3_states) { /* avoid accidental clash */
				if (a->id == 'V') {
					chi2 = chi1 + 2.0*M_PI_3;
				} 
				else {
					chi2 = chi1 - 2.0*M_PI_3;
				}
				adjust_to_closest_state(&chi2,a->id); /* rewrap within 2 PI */
			}
			a->chi2 = chi2;
		}
		else if (a->etc & G2_) {
			sprintf(error_string,"Amino acid %c%d has a 2nd gamma atom, despite not being ILE, THR or VAL!",a->id,a->num);
			stop(error_string);
		}
		else a->chi2 = chi2;
	} 
	else { /* check that G,A,P do not have gamma atoms */
		if (a->etc & G__ || a->etc & G2_) {
			sprintf(error_string,"Amino acid %c%d has gamma atom(s), despite being GLY, ALA or PRO!",a->id,a->num);
			stop(error_string);
		}
		a->chi2 = chi2;
		a->chi1 = chi1;
	}
}

/* convert Ramachandran angles to the relative triplet orientation */
void ramaset(triplet x, double phi, double psi, double tau, int cis)
{
	double *nn, sk;
	matrix t, tt, ttt;
	vector a = { 1.0, 0.0, 0.0 };

	if (cis) {
		sk = skew_;
		nn = n1_;
	} else {
		sk = skew;
		nn = n1;
		phi += M_PI;	/* trans-peptide magic */
	}

	rotmatrix(tt, a, tau - sk);	/* valence angle */
	rotmatrix(t, c1, -psi);
	rotation(ttt, tt, t);
	rotmatrix(t, nn, -phi);
	rotation(x, t, ttt);
}

/* omega rotations are not used as a rule */
void omegaset(triplet x, double omega, int cis)
{
	double *a;

	if (cis)
		a = cn1_;
	else {
		if (omega > 0)
			omega -= M_PI;
		else
			omega += M_PI;
		a = cn1;
	}

	rotmatrix(x, a, -omega);
}

/* find the best peptide bond orientation between two amino acids */
void amidorient(triplet x, AA *a, AA *b)
{
	double alpha, beta, gamma;
	vector zz, pp;

	if (a == NULL)		/* N-terminus */
		subtract(zz, b->ca, b->n);
	else if (b == NULL) {	/* C-terminus */
		subtract(zz, a->c, a->ca);
		subtract(pp, a->c, a->o);
		lincomb(zz, c[1] - o[1], zz, -c[1], pp);	/* inferring */
	} else if (b->etc & CIS) {	/* cis-peptide */
		subtract(zz, b->ca, a->ca);
		subtract(pp, b->n, a->ca);
		lincomb(zz, -ca_[1] - n_[1], zz, ca_[1], pp);	/* inferring */
	} else {		/* usual trans-peptide */
		subtract(zz, b->ca, a->ca);
		subtract(pp, b->n, a->c);
	}

	alpha = atan2(zz[0], -zz[1]);
	beta = atan2(sqrt(zz[0] * zz[0] + zz[1] * zz[1]), zz[2]);

	if (a == NULL) {
		beta -= atan2(n[1], -n[2]);	/* inferring */
		gamma = 0.0;
	} else {
		double zz2, pz;

		zz2 = square(zz);
		pz = dotprod(pp, zz);

		gamma = atan2(sqrt(zz2) * (pp[0] * zz[1] - pp[1] * zz[0]),
			      zz2 * pp[2] - pz * zz[2]);
	}

	eulerset(x, alpha, beta, gamma);
}

void aat_init(Chain * chain, Chaint * chaint){
  if(sizeof(chaint)->aat != chain->NAA * sizeof(AA)){
    (chaint)->aat = (AA *) realloc((chaint)->aat, chain->NAA * sizeof(AA));
    (chaint)->xaat = (triplet *) realloc((chaint)->xaat, chain->NAA * sizeof(triplet));
    (chaint)->ergt = (double *) realloc((chaint)->ergt, 5 * chain->NAA * sizeof(double));
    int i;	
    for (i = 1; i < chain->NAA; i++) {
	  (chaint)->aat[i].id = chain->aa[i].id;
	  (chaint)->aat[i].etc = chain->aa[i].etc;
	  (chaint)->aat[i].num = chain->aa[i].num;
	  (chaint)->aat[i].chainid = chain->aa[i].chainid;
	  (chaint)->aat[i].SCRot = chain->aa[i].SCRot;
    }
  }	
  if(sizeof(chaint)->xaat_prev != (chain->Nchains+1) * sizeof(triplet)){
    (chaint)->xaat_prev = (triplet *) realloc((chaint)->xaat_prev, (chain->Nchains+1) * sizeof(triplet));
  }
}

/* given the Ca coordinates and the amid orientation vectors {xaa(i)} (and xaa_prev for chain starts)
   (and perhaps the side chain dihedral angles), build the peptide */
void fulfill(Chain * chain,Chaint* chaint, simulation_params *sim_params)
{
	int i;

//fprintf(stdout,"FULFILL CALLED\n");
	/* calculate all backbone and side chain coordinates */
	for (i = 1; i < chain->NAA; i++)
		if (chain->aa[i].chainid != chain->aa[i-1].chainid) {
			acidate((chain->aa) + i, chain->xaa_prev[chain->aa[i].chainid], chain->xaa[i], sim_params);
		} else {
			acidate((chain->aa) + i, chain->xaa[i - 1], chain->xaa[i], sim_params);
		}

	chain->aa[0].id = -1;		/* dummy assignment at virtual position */
   
	aat_init(chain, chaint);
	
    
}

/* this allocates the memory for the data structures */
void allocmem_chain(Chain *chain, int NAA, int Nchains)
{
//	fprintf(stderr,"allocating chain for NAA=%d amino acids\n",NAA);
	(chain)->NAA = NAA;
	(chain)->Nchains = Nchains;
	(chain)->aa =  (AA*)realloc((chain)->aa, (chain)->NAA * sizeof(AA));
//	fprintf(stderr,"allocating erg: %ld\n",(chain)->NAA * (chain)->NAA * sizeof(double));
	(chain)->erg = (double*)realloc((chain)->erg, (chain)->NAA * (chain)->NAA * sizeof(double));
	///fprintf(stderr," %g",chain->erg);
	(chain)->xaa = (triplet*)realloc((chain)->xaa,(chain)->NAA * sizeof(triplet));
	(chain)->xaa_prev = (triplet*)realloc((chain)->xaa_prev,((chain)->Nchains + 1) * sizeof(triplet));
	if ((chain)->erg == NULL || (chain)->xaa == NULL || (chain)->aa == NULL || (chain)->xaa_prev == NULL) {
		if ((chain)->xaa == NULL) stop("allocmem_chain: Insufficient memory (chain->xaa)");
		if ((chain)->aa == NULL) stop("allocmem_chain: Insufficient memory (chain->aa)");
		if ((chain)->xaa_prev == NULL) stop("allocmem_chain: Insufficient memory (chain->xaa_prev)");
		if ((chain)->erg == NULL) stop("allocmem_chain: Insufficient memory (chain->erg)");
	}
}


/* initial configuration based on the provided sequence:
alpha-helix is in upper case, extended coil (beta-strand) is in lower case */
/* underscore means chain break (this is safe using utf-8 + various codes for helix/strand/etc), but building it as if it were a continuous chain */
void build_peptide_from_sequence(Chain * chain, Chaint *chaint, char *str, simulation_params *sim_params)
{
	int i;
	triplet helix, strand;

	matrix t = { { 1., 0., 0. },{ 0., 1., 0. },{ 0., 0., 1. } };

	int randx = 0, randy = 0, randz = 0;
	//if external AD grid choose the 8 corners of the box as starting position
	if ((sim_params->protein_model).external_potential_type == 5) {
		//fprintf(stderr, "random number %g %g %g \n", (double)rand(), RAND_MAX, sim_params->seed);
		randx = (round((double)rand() / RAND_MAX) * 2 - 1);
		randy = (round((double)rand() / RAND_MAX) * 2 - 1);
		randz = (round((double)rand() / RAND_MAX) * 2 - 1);
	}
	/* let's get rid of the separator, '_'
	   and also take notes of which chain the amino acid is in */
	fprintf(stderr,"Building protein from sequence: %s.\n",str);
	copy_string(&(sim_params->sequence),str);
	char * str_without_separator = (char *)malloc((strlen(str)+1)*sizeof(char));
	int * chain_ids = (int *)malloc((strlen(str)+2)*sizeof(int));
	int next = 0;
	int next_chain_id = 1; /* starting from 1 */
	for (int j = 0; str[j] != '\0'; j++){
		if (str[j]!='_') { /* copy amino acid characters */
			str_without_separator[next] = str[j]; //starting from 0
			chain_ids[next] = next_chain_id; //starting from 1
			next ++;
		} else { /* separator means new chain */
			next_chain_id ++;
		}
		//fprintf(stderr,"%c %d\n",str[j],next_chain_id);
	}
	str_without_separator[next] = '\0';
	str_without_separator = (char *) realloc(str_without_separator, (strlen(str_without_separator)+1)*sizeof(char));
	chain_ids = (int *) realloc(chain_ids, (strlen(str_without_separator))*sizeof(int));
	fprintf(stderr,"Building protein from sequence: %s (%d amino acids, %d chains).\n",str_without_separator, (int)strlen(str_without_separator), chain_ids[(int)strlen(str_without_separator)-1]);

	//fprintf(stderr,"Building peptide from sequence. ");
	//NAA = strlen(str) + 1;
	allocmem_chain(chain,strlen(str_without_separator)+1,next_chain_id);
	/* setting default parameters for the 0th amino acid (not used) */
	chain->aa[0].id = 'A';
	chain->aa[0].num = 0;
	chain->aa[0].chainid=0;
	chain->aa[0].etc=0;

	//sim_params->NAA = NAA;
	chain->Nchains = next_chain_id;
	chain->xaa_prev = malloc((next_chain_id+1)*sizeof(triplet));
      

	/* parse the string */
	for (i = 0; str_without_separator[i] != '\0'; i++) {
		chain->aa[i + 1].id = (str_without_separator[i] & COD) | 0x40;
		chain->aa[i + 1].etc = (str_without_separator[i] & ~COD);
		chain->aa[i + 1].num = i + 1;	/* enumerate */
		chain->aa[i + 1].chainid = chain_ids[i]; /* starting from 0 */
		fputc(str_without_separator[i], stderr);
	}
	if (chain_ids) free(chain_ids);

	/* Pauling and Corey alpha-helix with phi = -57 and psi = -45 */
	eulerset(helix, -1.25, 1.65, 2.15);

	/* right-twisted beta-strand with phi = -123 and psi = 135 */
	eulerset(strand, -1.7, 1.0, -1.2);

	transset(chain->xaa[0], t);
	transset(chain->xaa_prev[1], t);


	for (i = 1; i < chain->NAA; i++) {
		//chain->aa[i].SCRot = 0;
		

		if (chain->aa[i].etc & PSI)
			rotation(chain->xaa[i], t, strand);
		else
			rotation(chain->xaa[i], t, helix);

		//monitor chain breaks
		if (chain->aa[i].chainid != chain->aa[i-1].chainid) {
			//set xaa_prev for new chain
			casttriplet(chain->xaa_prev[chain->aa[i].chainid], chain->xaa[i]);
		}

		transset(t, chain->xaa[i]);
		/* fix accumulating errors every 64 residues */
		if ((i & 0x3F) == 0)
			fixtriplet(t);
	}

	//fprintf(stderr,"XAA\n");
	//for (i = 0; i < chain->NAA; i++) {
	//	fprintf(stderr,"%d ",i);
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");
	//fprintf(stderr,"XAA_PREV\n");
	//for (i = 1; i <= chain->Nchains; i++) {
	//	fprintf(stderr,"%d ",i);
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain->xaa_prev[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA_PREV\n");
	chain->aa[0].ca[0] = chain->aa[0].ca[1] = chain->aa[0].ca[2] = 0.0;
	
	if ((sim_params->protein_model).external_potential_type == 5 && transPtsCount!=0) {
		//srand(sim_params->seed);
		int transPtsID = rand() % transPtsCount;
		chain->aa[0].ca[0] = Xpts[transPtsID];
		chain->aa[0].ca[1] = Ypts[transPtsID];
		chain->aa[0].ca[2] = Zpts[transPtsID];
		//chain->aa[0].ca[0] = centerX - randx * (NX - 1) * spacing / 4;
		//chain->aa[0].ca[1] = centerY - randy * (NY - 1) * spacing / 4;
		//chain->aa[0].ca[2] = centerZ - randz * (NZ - 1) * spacing / 4;
		//fprintf(stderr, "random number %g %g %g \n", (double)rand(), RAND_MAX, sim_params->seed);
		//chain->aa[0].ca[0] = chain->aa[0].ca[1] = chain->aa[0].ca[2] = 0.0;
	}
	else {
		chain->aa[0].ca[0] = chain->aa[0].ca[1] = chain->aa[0].ca[2] = 0.0;	/* arbitrary */
	
	}




	for (i = 1; i < chain->NAA; i++) {
		//at the moment, xaa_prev is the same as xaa[i-1], so let's just use that
		carbonate_f((chain->aa) + i, (chain->aa) + i - 1, chain->xaa[i - 1]);
	}

	/* input diagnostics: sequence and length */
	fprintf(stderr, " %d\n", chain->NAA - 1);

	/* possibly initialise the sidechain dihedral angles.  It has to be done now,
	   because those cannot be calculated from only the backbone coordinates,
	   which we have when doing a crankshaft move */
	if ((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) {
	    for (i = 1; i < chain->NAA; i++) {
		initialise_sidechain_dihedral_angles((chain->aa)+i, &(sim_params->protein_model));
	    }
	}

	/* set existing atoms in aa.etc */
	/* this will be needed for chaint as well */
	for (int i=1; i<chain->NAA; i++) {
	    if (chain->aa[i].id == 'G') {
	        chain->aa[i].etc |= N__ | CA_ | C__ | O__ | H__ ;
	    } else if (sim_params->protein_model.use_gamma_atoms != NO_GAMMA) {
		/* gamma atom(s) */
		if (chain->aa[i].id == 'I' || chain->aa[i].id == 'T' || chain->aa[i].id == 'V') {
		    chain->aa[i].etc |= N__ | CA_ | C__ | O__ | CB_ | G__ | G2_ | H__ ;
		} else {
		    if (chain->aa[i].id == 'A') {
			chain->aa[i].etc |= N__ | CA_ | C__ | O__ | CB_ | H__ ;
		    } else {
			if (chain->aa[i].id == 'P') {
#ifdef LINUS_1995
			    /* LINUS 1995 doesn't have CG for PRO */
			    chain->aa[i].etc |= N__ | CA_ | C__ | O__ | CB_ ;
#else
			    chain->aa[i].etc |= N__ | CA_ | C__ | O__ | CB_ | G__ ;
#endif
			} else {
			    chain->aa[i].etc |= N__ | CA_ | C__ | O__ | CB_ | G__ | H__ ;
			}
		    }
		}
	    } else {
		chain->aa[i].etc |= N__ | CA_ | C__ | O__ | CB_ | H__ ;
	    }
	}

	/* build side chain and initialise energy */
	fulfill(chain,chaint,sim_params);

	free(str_without_separator);
}


/***********************************************************/
/****                MORE GEOMETRY TESTS                ****/
/***********************************************************/


/* Check if amino acid a is an h-bond donor for amino acid b,
   by comparing the d(O,H) distance and the <(O,H,N) and <(C,O,H) angles
   to the cutoff values. */
int hdonor(AA *a, AA *b, model_params *mod_params)
{
	vector oh;

	if (!((a->etc & H__) && (b->etc & O__))) return 0;

	subtract(oh, a->h, b->o);
	if (square(oh) < mod_params->hboh2) {
		vector hn;
		subtract(hn, a->n, a->h);
		if (cosine(oh, hn) > mod_params->hbohn) {
			vector co;
			subtract(co, b->o, b->c);
			if (cosine(co, oh) > mod_params->hbcoh) {
				return 1;
			}
		}
	}
	return 0;
}

/* Calculate the strength of hbond, using distance and angle hard cutoff functions.
   This is the number, which the H bond strength constant is multiplied with. 
*/
/* double hstrength_hard_cutoff(vector N, vector H, vector O, vector C, model_params *mod_params){
        double ans = 1.0;
        vector oh;
        subtract(oh, H, O);
        vector hn;
        subtract(hn, N, H);
        vector co;
        subtract(co, O, C);
        double cos_oh_hn = cosine(oh,hn);
        if(cos_oh_hn < mod_params->hbohn) return 0;
        double cos_co_oh = cosine(co,oh);
        if(cos_co_oh < mod_params->hbcoh) return 0;
        double oh_2 = square(oh);
        if(oh_2 > mod_params->hboh2) return 0; 
        return ans;
} */

/* Calculate the strength of hbond, using distance and angle linear cutoff functions.
   This is the number, which the H bond strength constant is multiplied with.
*/
/* inline double linear_decay_local(double distance,
		           double contact_cutoff,
			   //double strength,
			   double decay_width) {

	if (distance > contact_cutoff + decay_width) return 0.0;
	if (distance < contact_cutoff) return 1.0;
	return (distance - contact_cutoff) / decay_width;
}
double hstrength_linear(vector N, vector H, vector O, vector C, model_params *mod_params){
	double ans = 1.0;
	vector oh;
	subtract(oh, H, O);
	vector hn;
	subtract(hn, N, H);
	vector co;
	subtract(co, O, C);
	double cos_oh_hn = cosine(oh,hn);
	if(cos_oh_hn < mod_params->hbohn-mod_params->hbohn_decay_width) return 0;
	double cos_co_oh = cosine(co,oh);
	if(cos_co_oh < mod_params->hbcoh-mod_params->hbcoh_decay_width) return 0;
	double oh_2 = square(oh);
	if(oh_2 > mod_params->hboh2) {
		double abs_oh = sqrt(oh_2);
		double abs_cut = sqrt(mod_params->hboh2);
		ans *= linear_decay_local(abs_oh, abs_cut, mod_params->hboh_decay_width);
//fprintf(stderr,"doh %g %g %g\n",abs_oh,abs_cut,linear_decay_local(abs_oh, abs_cut, mod_params->hboh_decay_width));
	}
	if(cos_oh_hn < mod_params->hbohn) ans *= (mod_params->hbohn - cos_oh_hn) / mod_params->hbohn_decay_width;
	if(cos_co_oh < mod_params->hbcoh) ans *= (mod_params->hbcoh - cos_co_oh) / mod_params->hbcoh_decay_width;
//fprintf(stderr,"ohn %g %g %g\n",cos_oh_hn,mod_params->hbohn,(mod_params->hbohn - cos_oh_hn) / mod_params->hbohn_decay_width);
//fprintf(stderr,"coh %g %g %g\n",cos_oh_hn,mod_params->hbcoh,(mod_params->hbcoh - cos_co_oh) / mod_params->hbcoh_decay_width);
		
	return ans;
} */


/* Calculate the strength of hbond, using distance and angle cutoff functions.
   This is the number, which the H bond strength constant is multiplied with. 
   The cutoff functions are:
       ( d_cut(O,H) / d(O,H) )^4,           if d(O,H) > d_cut(O,H)
       ( cos (OH,HN) / cos_cut(OH,HN) )^4,    if cos(OH,HN) < cos_cut(OH,HN)
       ( cos (CO,OH) / cos_cut(CO,OH) )^4,    if cos(CO,OH) < cos_cut(CO,OH)
*/
double hstrength(vector N, vector H, vector O, vector C, model_params *mod_params){
	double ans = 1.0;
	vector oh;
	subtract(oh, H, O);
	vector hn;
	subtract(hn, N, H);
	vector co;
	subtract(co, O, C);
	double cos_oh_hn = cosine(oh,hn);
	if(cos_oh_hn < 0) return 0;
	double cos_co_oh = cosine(co,oh);
	if(cos_co_oh < 0) return 0;
	if(square(oh) > mod_params->hboh2) ans *= (mod_params->hboh2/square(oh))*(mod_params->hboh2/square(oh));
	if(cos_oh_hn < mod_params->hbohn) ans *= (cos_oh_hn/mod_params->hbohn)*(cos_oh_hn/mod_params->hbohn)*(cos_oh_hn/mod_params->hbohn)*(cos_oh_hn/mod_params->hbohn) ;
	if(cos_co_oh < mod_params->hbcoh) ans *= (cos_co_oh/mod_params->hbcoh)*(cos_co_oh/mod_params->hbcoh)*(cos_co_oh/mod_params->hbcoh)*(cos_co_oh/mod_params->hbcoh) ;
		
	return ans;
}

/* Check if amino acids a and b are in contact,
   that is the CA or CB atoms of d(a,b) < touch2.
   This is a very crude algorithm, and it does not work for gamma atoms. */
int contact(AA *a, AA *b, model_params *mod_params)
{
	vector x, y;

	lincomb(x, 1. - mod_params->part, a->ca, mod_params->part, a->cb);
	lincomb(y, 1. - mod_params->part, b->ca, mod_params->part, b->cb);

	if (distance(x, y) > mod_params->touch2)
		return 0;

	if (-1.0 < mod_params->split && mod_params->split < 1.0) {
		subtract(x, a->cb, a->ca);
		subtract(y, b->cb, b->ca);

		if (cosine(x, y) < mod_params->split)
			return 0;
	}

	return 1;
}

/* Check if the side chains of two amino acids point in the same direction. */
int aligned( AA *a, AA *b)
{
	vector x, y;

	subtract(x, a->cb, a->ca);
	subtract(y, b->cb, b->ca);

	return cosine(x, y) > -0.11111111111111111111;
}

/* calculate phi and psi angles relying on triplets instead of neighbors */
/*static double phi_angle(struct AA *a, triplet xnt)
{
	vector nca, cac;	// N-Ca and Ca-C bonds 

	subtract(nca, a->ca, a->n);
	subtract(cac, a->c, a->ca);

	return dihedral(xnt[1], nca, cac);
}

static double psi_angle(struct AA *a, triplet xct)
{
	vector nca, cac;	// N-Ca and Ca-C bonds 

	subtract(nca, a->ca, a->n);
	subtract(cac, a->c, a->ca);

	return dihedral(nca, cac, xct[1]);
}*/


/***********************************************************/
/****            AMINO ACID AND PEPTIDE  I/O            ****/
/***********************************************************/

/****                     PRINTING                      ****/

/* Printing an amino acid residue in PDB format */
int pdbrecord( AA *a, int j, model_params *mod_params, FILE *outfile)
{
	char fmt[] = "ATOM  %5d %4s XAA A9999    %8.3f%8.3f%8.3f\n";
	//const char chain = 'A';
	char *gatom = "testingtesting";
	char *g2atom = "testingtesting";

	//fprintf(stderr,"%d ",a->chainid);
	sprintf(fmt + 14, "%3s %c%4d", aa123(a->id), 'A' + a->chainid - 1, a->num & 0xFFF);
	fmt[23] = ' ';

	fprintf(outfile,fmt, ++j, " N  ", a->n[0], a->n[1], a->n[2]);
	fprintf(outfile,fmt, ++j, " CA ", a->ca[0], a->ca[1], a->ca[2]);
	fprintf(outfile,fmt, ++j, " C  ", a->c[0], a->c[1], a->c[2]);
	fprintf(outfile,fmt, ++j, " O  ", a->o[0], a->o[1], a->o[2]);
	if (a->id != 'G') {
		fprintf(outfile,fmt, ++j, " CB ", a->cb[0], a->cb[1], a->cb[2]);
		if (mod_params->use_gamma_atoms != NO_GAMMA) {
#ifdef LINUS_1995
			/* LINUS 1995 doesn't have CG for PRO, so we have to skip for PRO */
			if (a->id != 'A' && a->id != 'P' && a->g[2] < FAR) {
#else
			if (a->id != 'A' && a->g[2] < FAR) {
#endif
				switch (a->id) {
				case 'I':
				case 'V':
					gatom = " CG1";
					g2atom = " CG2";
					break;
				case 'T':
					gatom = " OG1";
					g2atom = " CG2";
					break;
				case 'S':
					gatom = " OG ";
					break;
				case 'C':
					gatom = " SG ";
					break;
				case 'U':
					gatom = "SEG ";
					break;
				case 'X':
					/* TODO: would it ever get here? */
					fprintf(stderr,"It should not ever get here.");
					break;
				default:
					gatom = " CG ";
				}
				if (a->etc & G__) fprintf(outfile,fmt, ++j, gatom, a->g[0], a->g[1], a->g[2]);
				if ((a->id == 'I' || a->id == 'V' || a->id == 'T') && a->g2[2] < FAR) {
					if (a->etc & G2_) fprintf(outfile,fmt, ++j, g2atom, a->g2[0], a->g2[1], a->g2[2]);
				}
			}
		}
	}

	if (a->id != 'P')
		fprintf(outfile,fmt, ++j, " H  ", a->h[0], a->h[1], a->h[2]);

	return j;
}

/* Printing an amino acid chain in PDB format */
void pdbprint( AA *a, int count, model_params *mod_params, FILE *outfile, double *E_tot)
{
	int i, j;
	static int model = 1;

	fprintf(outfile,"MODEL     %4d\n", model);
	if (E_tot != NULL) fprintf(outfile,"REMARK ENERGY %.6f\n", *E_tot);

	for (j = 0, i = 1; i < count; i++) {
		if (i>1 && a[i].chainid != a[i-1].chainid) {
			fprintf(outfile,"TER   %5d      %3s %c%4d\n",
			       j, aa123(a[i].id), 'A', i);
		}
		j = pdbrecord(a + i, j, mod_params, outfile);
	}

	fprintf(outfile,"TER   %5d      %3s %c%4d\n",
	       ++j, aa123(a[i - 1].id), 'A', i - 1);
	fprintf(outfile,"ENDMDL\n");
	model++;
}

/****                      READING                      ****/

/* Increase the size of allocated memory for the amino acid array pointer */
static int dblalloc(AA **paa, int *size)
{

  
	AA *tmp = NULL;
	int ns = 2 * (*size + 1);

   
    tmp = ( AA *) realloc(*paa, ns * sizeof(AA));
       
    

	if (tmp == NULL)
		return 0;

	*size = ns;
	*paa = tmp;
	
	
	
	
	return ns;
}

/* Read an amino acid residue from a PDB input, possibly with missing atoms.
   Reading it line-by-line; new residue number.  ENDMDL causes interruption.
   Works for multi-chain proteins.  TER, jump in the residue number
   or change in the chain id cause the chainid counter to advance. */
int getaa(AA *a, FILE *infile)
{
	static char line[83] = "Veni Vidi Vici";
	int mask, status = 0, pos = 0x7FFF, cur;
	double *atm;
	char prev_chain = ' ';
	int k;

	a->etc = 0x0;
	a->etc &= ~ATM;		/* clean up */
	a->id = 'X';

	do {
		//fprintf(stderr,"GETAA: %s\n",line);
		/* break at ENDMDL or END */
		if (line[0] == 'E' && line[3] == 'M' && line[4] == 'D') {
			//advance the file pointer to the next line
//			fprintf(stderr,"ENDMDL reached %d. return\n",(int)0x7FFF);
			fgets(line, sizeof(line), infile);
			//fprintf(stderr,"GETAA cheat: %s\n",line);
			pos = 0x7FFF;
			return pos; //TODO check it is OK to return
		}
		if (line[0] == 'E' && line[1] == 'N' && line[2] == 'D') {
//			fprintf(stderr,"END reached %d. return\n",0x7FFF);
			fgets(line, sizeof(line), infile);
			//fprintf(stderr,"GETAA cheat: %s\n",line);
			pos = 0x7FFF;
			return pos; //TODO check it is OK to return
		}

		/* chain break at TER */
		if (line[0] == 'T' && line[1] == 'E' && line[2] == 'R') {
//			fprintf(stderr,"TER reached %d. not return\n",(int)0x7FFE);
			pos = 0x7FFE; //chain break
		}

		/* check for ATOM and HETATM record */
		if ((line[0] != 'A' || line[1] != 'T' || line[3] != 'M') &&
		    (line[0] != 'H' || line[3] != 'A' || line[4] != 'T')) {
			//fprintf(stderr,"line to ignore: %s\n",line);
			continue;
		}

		/* check for change in the chain name (22) */
//fprintf(stderr,"%c",line[21]);
		if (status == 0) prev_chain = line[21];
		if (status && line[21] != prev_chain) {//chain break
			fprintf(stderr,"Found a chain break (%c->%c) return\n",prev_chain,line[21]);
			return 0x7FFE;
		}

		/* check for change in the residue number (23:26) up to 4 digits */
		if ((k = sscanf(line+22,"%4d",&cur)) != 1){
			stop("Could not read the amino acid number\n");
		//} else {
		//	fprintf(stderr,"Amino acid %d\n",cur);
		}
		if (status && cur != pos) {
			if (cur - pos == 1) {
				//next amino acid
				//fprintf(stderr,"next amino acid return\n");
				return pos; //!= 0x7FFF;
			} else {
				//chain gap, i.e. end of chain
				fprintf(stderr,"Found a chain break (chain %c %d->%d) return\n",prev_chain,pos,cur);
				return 0x7FFE;
			}
		}
		pos = cur;
		status = 1;

		/* check the atom type */
		if (line[13] == 'N' && line[14] == ' ') {
			atm = a->n;
			mask = N__;
		} else if (line[13] == 'C' && line[14] == 'A') {
			atm = a->ca;
			mask = CA_;
		} else if (line[13] == 'C' && line[14] == ' ') {
			atm = a->c;
			mask = C__;
		} else if (line[13] == 'O' && line[14] == ' ') {
			atm = a->o;
			mask = O__;
		} else if (line[13] == 'C' && line[14] == 'B') {
			atm = a->cb;
			mask = CB_;
		} else if (line[13] == 'H' && line[14] == ' ') {
			atm = a->h;
			mask = H__;
		} else if (line[13] != 'H' && line[14] == 'G'
			   && line[15] != '2') {
#ifdef LINUS_1995
			/* LINUS 1995 doesn't have CG for PRO */
			if (line[17]=='P' && line[18]=='R' && line[19]=='O') {
				//fprintf(stderr, "Skipping proline CG\n");
				continue;
			}
#endif
			atm = a->g;
			mask = G__;
		} else if (line[13] != 'H' && line[14] == 'G'
			   && line[15] == '2') {
			atm = a->g2;
			mask = G2_;
		} else
			continue;

		/* keep the first alternate location */
		if ((mask & a->etc) && line[17] != ' ')
			continue;

		/* scan coordinates and check for a good scan */
		if (sscanf(line + 30, "%8lf%8lf%8lf", atm, atm + 1, atm + 2)
		    != 3)
			continue;

		if (CA_ & mask & ~a->etc) {
			a->num = cur;	/* may have duplicates */
			if ((a->id = aa321(line + 17)) == 'X')
				fprintf(stderr, "Strange residue  %.9s\n",
					line + 17);
		}

		/* keep track of what's read */
		a->etc |= mask;
	} while (fgets(line, sizeof(line), infile) != NULL);

	fprintf(stderr,"eof return\n");
	return EOF;
}

/* read a PDB input residue after residue, store coordinates in an array of
structures, return the number of successfully read amino acids, possibly with
missing atoms. At most one less than size amino acids are completely read in.
The storage is doubled on reaching its size. */
int getpdb(AA **ptraa, int *size, int *Nchains, FILE *infile)
{
	
	int pos, status, tot = 1;
	AA b;
	int next_chainid = 1;
   
	do {
		
		/* keep track of memory */
		if (tot >= *size && !dblalloc(ptraa, size)) {
			stop("getpdb: Insufficient memory in getpdb");
		}
       
		pos = getaa(&b, infile);
		//fprintf(stderr,"pos status: %d\n",pos);

		if (b.etc & CA_) {	/* advance on CA only */
			/* mark chain id */
			b.chainid = next_chainid;
			(*ptraa)[tot] = b;
			tot++;
		}

		/* set status from pos */
		if (pos == EOF) {
//			fprintf(stderr,"End of file.\n");
			status = 0;
		} else if (pos == 0x7FFF) { /* ENDMDL */
			/* end of molecule, stop */
//			fprintf(stderr,"End of molecule.\n");
			status = 0;
		} else {
			/* reset status, keep going */
			status = 1;
			if (pos == 0x7FFE) { /* TER or chain break */
//				fprintf(stderr,"End of chain.\n");
				next_chainid ++;
			}
		}
		
	} while (status == 1);
	if (tot > 1)
		fprintf(stderr,"getpdb: read in %d chains (NAA=%d)\n",next_chainid,tot-1);
	*Nchains = next_chainid;

	if (tot > 1) {
		*size = tot;
		return tot - 1;
	} else
		return EOF;
}

/***********************************************************/
/****              SETUP  / INPUT / OUTPUT              ****/
/***********************************************************/


/* this frees the memory for the data structures */
void freemem_chain(Chain *chain)
{
	if(chain){
	    if (chain->aa ) {
		free(chain->aa);
		chain->aa = NULL;
	    }
	    if (chain->erg) {
		free(chain->erg);
		chain->erg = NULL;
	    }
	    if (chain->xaa) {
		free(chain->xaa);
		chain->xaa = NULL;
	    }
	    if (chain->xaa_prev) {
		free(chain->xaa_prev);
		chain->xaa_prev = NULL;
	    }
	}
}

void freemem_chaint(Chaint *chaint)
{
	if(chaint){
	    if (chaint->aat) {
		free(chaint->aat);
		chaint->aat = NULL;
	    }
	    if (chaint->ergt) {
		free(chaint->ergt);
		chaint->ergt = NULL;
	    }
	    if (chaint->xaat) {
		free(chaint->xaat);
		chaint->xaat = NULL;
	    }
	    if (chaint->xaat_prev) {
		free(chaint->xaat_prev);
		chaint->xaat_prev = NULL;
	    }
    } 
}

/* minimally adjust coordinates to better equalize Ca-Ca distances */
static double repair_segment( Chain *chain1, Chain *chain2, int aa_start, int aa_end)
{
	int i;
	double p, rmse = 0.0;
	extern const double caca, caca_;
	double *diag, *cosn, *dl, *ilen;
	vector a;

	int NAA = aa_end - aa_start + 2;
	int chainid = chain1->aa[aa_start].chainid;

	fprintf(stderr,"REPAIR SEGMENT %d -- %d\n",aa_start,aa_end);

	//fprintf(stderr,"XAA\n");
	//for (i = 1; i < NAA; i++) {
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain1->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");

	/* borrow some storage */
	diag = (double *) malloc(sizeof(double) * NAA * 4);
	cosn = diag + NAA;
	dl = diag + 2 * NAA;
	ilen = diag + 3 * NAA;
	for (i = 0; i < NAA * 4 ; i++) {
		diag[i] = 0.0;
	}

	/* evaluate directional unit vectors and adjustments */
	for (i = 1; i < NAA - 1; i++) {
		subtract(a, chain2->aa[aa_start-1 + i + 1].ca, chain2->aa[aa_start-1+i].ca);
		p = sqrt(square(a));

		dl[i] = ((chain2->aa[aa_start-1+i + 1].etc & CIS) ? caca_ : caca) - p;

		p = 1. / p;
		scale(chain1->xaa[aa_start-1+i][2], p, a);	/* normalize */

		ilen[i] = p;	/* save for later */
		diag[i] = -2.0;	/* main diagonal */
	}

	//fprintf(stderr, "diag:\n");
	//for (i = 0; i < NAA * 4 ; i++) {
	//	fprintf(stderr, "%g ",diag[i]);
	//	if (i%chain1->NAA == 0){
	//		fprintf(stderr, "\n");
	//	}
	//}
	//fprintf(stderr, "\n");

	//fprintf(stderr,"XAA\n");
	//for (i = 1; i < NAA; i++) {
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain1->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");


	/* evaluate off-diagonal directional cosines */
	for (i = 1; i < NAA - 2; i++)
		cosn[i] = dotprod(chain1->xaa[aa_start-1+i][2], chain1->xaa[aa_start-1+i + 1][2]);

	/* evaluate Lagrange multipliers by Gaussian elimination */
	for (i = 1; i < NAA - 2; i++) {
		p = cosn[i] / diag[i];
		diag[i + 1] -= cosn[i] * p;
		dl[i + 1] -= dl[i] * p;
		cosn[i] = p;	/* reused in the next loop */
	}

	for (i = NAA - 3; i > 0; i--)
		dl[i] -= dl[i + 1] * cosn[i];

	for (i = 1; i < NAA - 1; i++)
		dl[i] /= diag[i];

	//fprintf(stderr, "diag:\n");
	//for (i = 0; i < NAA * 4 ; i++) {
	//	fprintf(stderr, "%g ",diag[i]);
	//	if (i%chain1->NAA == 0){
	//		fprintf(stderr, "\n");
	//	}
	//}
	//fprintf(stderr, "\n");

	/* adjust alpha-carbons */
	dl[0] = dl[NAA - 1] = 0.;
	chain1->xaa_prev[chainid][2][0] = chain1->xaa[aa_end][2][0] = 0.;
	chain1->xaa_prev[chainid][2][1] = chain1->xaa[aa_end][2][1] = 0.;
	chain1->xaa_prev[chainid][2][2] = chain1->xaa[aa_end][2][2] = 0.;
	for (i = 1; i < NAA; i++) {
		if (i==1) { //use xaa_prev for the start
			lincomb(chain1->xaa[aa_start-1+i][0], -dl[i - 1], chain1->xaa_prev[chainid][2], dl[i], chain1->xaa[aa_start-1+i][2]);
		} else {
			lincomb(chain1->xaa[aa_start-1+i][0], -dl[i - 1], chain1->xaa[aa_start-1+i - 1][2], dl[i], chain1->xaa[aa_start-1+i][2]);
		}
		rmse += square(chain1->xaa[aa_start-1+i][0]);

		add(chain1->aa[aa_start-1+i].ca, chain2->aa[aa_start-1+i].ca, chain1->xaa[aa_start-1+i][0]);
	}

	//fprintf(stderr,"XAA\n");
	//for (i = 1; i < NAA; i++) {
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain1->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");

	/* adjust nitrogens and carbons */
	add(chain1->aa[aa_start].n, chain2->aa[aa_start].n, chain1->xaa[aa_start][0]);
	for (i = 1; i < NAA - 1; i++) {
		subtract(a, chain2->aa[aa_start-1+i].c, chain2->aa[aa_start-1+i].ca);
		//fprintf(stderr,"a %g %g %g\n",a[0],a[1],a[2]);
		p = dotprod(a, chain1->xaa[aa_start-1+i][2]) * ilen[i];
		//fprintf(stderr,"p %g\n",p);
		lincomb(a, 1. - p, chain1->xaa[aa_start-1+i][0], p, chain1->xaa[aa_start-1+i + 1][0]);
		//fprintf(stderr,"a %g %g %g\n",a[0],a[1],a[2]);

		//fprintf(stderr,"chain2.c %g %g %g\n",chain2->aa[i].c[0],chain2->aa[i].c[1],chain2->aa[i].c[2]);
		add(chain1->aa[aa_start-1+i].c, chain2->aa[aa_start-1+i].c, a);
		//fprintf(stderr,"chain1.c %d %g %g %g\n",i,chain1->aa[i].c[0],chain1->aa[i].c[1],chain1->aa[i].c[2]);

		subtract(a, chain2->aa[aa_start-1+i + 1].ca, chain2->aa[aa_start-1+i + 1].n);
		p = dotprod(a, chain1->xaa[aa_start-1+i][2]) * ilen[i];
		//fprintf(stderr,"p %g\n",p);
		lincomb(a, p, chain1->xaa[aa_start-1+i][0], 1. - p, chain1->xaa[aa_start-1+i + 1][0]);
		//fprintf(stderr,"a %g %g %g\n",a[0],a[1],a[2]);

		//fprintf(stderr,"last chain1.n %d %g %g %g\n",i+1,chain1->aa[i+1].n[0],chain1->aa[i+1].n[1],chain1->aa[i+1].n[2]);
		add(chain1->aa[aa_start-1+i + 1].n, chain2->aa[aa_start-1+i + 1].n, a);
		//fprintf(stderr,"last chain1.n %d %g %g %g\n",i+1,chain1->aa[i+1].n[0],chain1->aa[i+1].n[1],chain1->aa[i+1].n[2]);
	}
	//fprintf(stderr,"? %d %d\n",NAA-1,aa_end);
	add(chain1->aa[aa_end].c, chain2->aa[aa_end].c, chain1->xaa[aa_end][0]);
	//fprintf(stderr,"last chain1.c %d %g %g %g\n",i,chain1->aa[i].c[0],chain1->aa[i].c[1],chain1->aa[i].c[2]);

	/* for (i = 1; i < NAA - 1; i++)
	   fprintf(stderr, "%f %f % f\n", len[i],
	   sqrt(distance(aa[i + 1].ca, aa[i].ca)), dl[i]);

	   fprintf(stderr, "RMSE = %g\n", sqrt(rmse / (NAA - 1))); */

	free(diag);

	return rmse;
}

/* minimally adjust coordinates to better equalize Ca-Ca distances */
// TODO: multi-chain proteins
static double repair_multichain( Chain *chain1, Chain *chain2) {

	int chain_starts[chain1->Nchains]; //numbering from 0
	int chain_ends[chain1->Nchains]; //numbering from 0

	//fprintf(stderr,"REPAIR MULTICHAIN\n");
	int next_chain = 0;
	chain_starts[next_chain] = 1;
	//fprintf(stderr,"start %d",1);
	for (int i = 2; i < chain1->NAA; i++) {
		if (chain1->aa[i].chainid != chain1->aa[i-1].chainid) {
			chain_ends[next_chain++] = i-1;
			chain_starts[next_chain] = i;
			//fprintf(stderr,", end %d, start %d",i-1,i);
		}
	}
	chain_ends[next_chain] = chain1->NAA-1;
	//fprintf(stderr,", end %d\n",chain1->NAA-1);
	if (next_chain+1 != chain1->Nchains) {
		fprintf(stderr,"next_chain = %d, Nchains = %d\n",next_chain+1,chain1->Nchains);
		stop("Something has gone wrong while finding chain starts and chain ends.\n");
	}

	double rmse = 0.0;
	for (int i = 0; i < chain1->Nchains; i++) {
		rmse += repair_segment(chain1, chain2, chain_starts[i], chain_ends[i]);
	}

	return rmse;

}

/* minimally adjust coordinates to better equalize Ca-Ca distances */
// TODO: multi-chain proteins
static double repair( Chain *chain1, Chain *chain2)
{
	int i;
	double p, rmse = 0.0;
	extern const double caca, caca_;
	double *diag, *cosn, *dl, *ilen;
	vector a;

	/* borrow some storage */
	diag = (double *) malloc(sizeof(double) * chain1->NAA * 4);
	cosn = diag + chain1->NAA;
	dl = diag + 2 * chain1->NAA;
	ilen = diag + 3 * chain1->NAA;

	/* mark chain starts and chain ends for multi-chain proteins,
	   we will adjust their CA only using 1 neighbour */
	int *chain_ends, *chain_starts;
	chain_starts = (int *) malloc(sizeof(int) * chain1->NAA);
	chain_ends = (int *) malloc(sizeof(int) * chain1->NAA);
	chain_starts[0] = 0; chain_ends[0] = 0;
	for (i = 1; i < chain1->NAA; i++) {
		if (i == 1) {
			chain_starts[i] = 1;
		} else {
			if (chain1->aa[i].chainid != chain1->aa[i-1].chainid)
				chain_starts[i] = 1;
			else
				chain_starts[i] = 0;
		}
		if (i == chain1->NAA -1) {
			chain_ends[i] = 1;
		} else {
			if (chain1->aa[i].chainid != chain1->aa[i+1].chainid)
				chain_ends[i] = 1;
			else
				chain_ends[i] = 0;
		}
	}

	/* evaluate directional unit vectors and adjustments */
	for (i = 1; i < chain1->NAA - 1; i++) {
		subtract(a, chain2->aa[i + 1].ca, chain2->aa[i].ca);
		p = sqrt(square(a));
		//CA-CA length difference from the ideal value
		dl[i] = ((chain2->aa[i + 1].etc & CIS) ? caca_ : caca) - p;

		//xaa[i][2] is the unit vector pointing from CA(i)->CA(i+1)
		p = 1. / p;
		scale(chain1->xaa[i][2], p, a);	/* normalize */

		//length of actual CA-CA bond
		ilen[i] = p;	/* save for later */
		diag[i] = -2.0;	/* main diagonal */
	}

	/* evaluate off-diagonal directional cosines */
	for (i = 1; i < chain1->NAA - 2; i++)
		cosn[i] = dotprod(chain1->xaa[i][2], chain1->xaa[i + 1][2]);

	/* evaluate Lagrange multipliers by Gaussian elimination */
	/* separately for all chains */
	for (i = 1; i < chain1->NAA - 2; i++) {
		//special cases
		if (chain_starts[i] == 1) //reset for start of the chain
			diag[i] = -2.0;
		if (chain_ends[i] != 1 && chain_ends[i+1] != 1) {//do not change the next parameters for chain ends
			p = cosn[i] / diag[i];
			diag[i + 1] -= cosn[i] * p;
			dl[i + 1] -= dl[i] * p;
			cosn[i] = p;	/* reused in the next loop */
		} else {
			fprintf(stderr,"Missing something at residue %d, near the chain break.\n",i);
		}
	}

//	fprintf(stderr,"dl ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",dl[i]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");
//	fprintf(stderr,"cosn ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",cosn[i]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");
//	fprintf(stderr,"diag ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",diag[i]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");
//	fprintf(stderr,"ilen ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",ilen[i]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");


	for (i = chain1->NAA - 3; i > 0; i--)
		if (chain_ends[i] != 1 && chain_ends[i+1] != 1) //do not change the next parameters for chain ends
			dl[i] -= dl[i + 1] * cosn[i];
		else
			fprintf(stderr,"Missing something at residue %d, too, near the chain break.\n",i);

	for (i = 1; i < chain1->NAA - 1; i++)
		dl[i] /= diag[i];

//	fprintf(stderr,"dl ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",dl[i]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");

	/* adjust alpha-carbons */
	dl[0] = dl[chain1->NAA - 1] = 0.;
	chain1->xaa[0][2][0] = chain1->xaa[chain1->NAA - 1][2][0] = 0.;
	chain1->xaa[0][2][1] = chain1->xaa[chain1->NAA - 1][2][1] = 0.;
	chain1->xaa[0][2][2] = chain1->xaa[chain1->NAA - 1][2][2] = 0.;
	for (i = 1; i < chain1->NAA; i++) {
		/* at chain ends only move using the neighbour on the same chain */
		double dl_own, dl_prev;
		if (chain_ends[i] == 1) //end of chain
			dl_own = 0;
		else
			dl_own = dl[i];
		if (chain_starts[i] == 1) //beginning of chain
			dl_prev = 0.;
		else
			dl_prev = -dl[i - 1];
		lincomb(chain1->xaa[i][0], dl_prev, chain1->xaa[i - 1][2], dl_own, chain1->xaa[i][2]);
		rmse += square(chain1->xaa[i][0]);

		add(chain1->aa[i].ca, chain2->aa[i].ca, chain1->xaa[i][0]);
	}
//	fprintf(stderr,"ca ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",chain1->aa[i].ca[0]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");
//	fprintf(stderr,"ca ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",chain1->aa[i].ca[1]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");
//	fprintf(stderr,"ca ");
//	for (i = 1; i < chain1->NAA; i++) {
//		fprintf(stderr,"%g ",chain1->aa[i].ca[2]);
//		if (chain_ends[i] == 1) fprintf(stderr,"\n");
//	}
//	fprintf(stderr,"\n");

	/* adjust nitrogens and carbons */
	for (i = 1; i < chain1->NAA; i++) {
		// amide N of its own amino acid for chain starts
		if (chain_starts[i] == 1) //for chain starts, just shift its own n with ca
			add(chain1->aa[i].n, chain2->aa[i].n, chain1->xaa[i][0]);
		// carbonil C
		if (chain_ends[i] == 1) { //for chain ends, just shift c with ca
			add(chain1->aa[i].c, chain2->aa[i].c, chain1->xaa[i][0]);
		} else {
			//carbonil C of the peptide bond
			subtract(a, chain2->aa[i].c, chain2->aa[i].ca);
			p = dotprod(a, chain1->xaa[i][2]) * ilen[i];
			lincomb(a, 1. - p, chain1->xaa[i][0], p, chain1->xaa[i + 1][0]);

			add(chain1->aa[i].c, chain2->aa[i].c, a);

			// amide N of the same peptide bond (N of the next amino acid)
			subtract(a, chain2->aa[i + 1].ca, chain2->aa[i + 1].n);
			p = dotprod(a, chain1->xaa[i][2]) * ilen[i];
			lincomb(a, p, chain1->xaa[i][0], 1. - p, chain1->xaa[i + 1][0]);

			add(chain1->aa[i + 1].n, chain2->aa[i + 1].n, a);
		}
	}

	/* for (i = 1; i < NAA - 1; i++)
	   fprintf(stderr, "%f %f % f\n", len[i],
	   sqrt(distance(aa[i + 1].ca, aa[i].ca)), dl[i]);

	   fprintf(stderr, "RMSE = %g\n", sqrt(rmse / (NAA - 1))); */

	free(diag);
	free(chain_starts);
	free(chain_ends);

	return rmse;
}

/* parsing ATOM entries from a PDB file for initialization
returning the number of parsed amino acids */
// this routine does not change the coordinates
int pdbin(Chain *chain, simulation_params *sim_params, FILE *infile)
{
	int retv;

//	fprintf(stderr,"Reading in PDB input from file.\n");

	/* scan the PDB-like input from infile, and allocate the memory */

    Chain *tempchain = malloc(sizeof(Chain));
	tempchain->NAA = 0; tempchain->Nchains = 0; tempchain->aa = NULL; tempchain->erg = NULL; tempchain->xaa = NULL; tempchain->xaa_prev = NULL;
	retv = getpdb(&(tempchain->aa), &(tempchain->NAA), &(tempchain->Nchains), infile);
//	fprintf(stderr,"tempchain->NAA=%d\n",tempchain->NAA);
	
	tempchain->xaa = (triplet*)realloc(tempchain->xaa, tempchain->NAA * sizeof(triplet));
	for (int i = 0; i < tempchain->NAA; i++) {
		for (int j=0; j<3; j++) {
			for (int k=0; k<3; k++) {
				tempchain->xaa[i][j][k] = 0.0;
			}
		}
	}
	tempchain->erg = (double*)realloc(tempchain->erg, tempchain->NAA * tempchain->NAA * sizeof(double));	
//	fprintf(stderr,"Allocating memory for xaa_prev (Nchains=%d)\n",tempchain->Nchains);
	tempchain->xaa_prev = (triplet*)realloc(tempchain->xaa_prev, (tempchain->Nchains+1) * sizeof(triplet));	
	for (int i = 0; i <= tempchain->Nchains; i++) {
		for (int j=0; j<3; j++) {
			for (int k=0; k<3; k++) {
				tempchain->xaa_prev[i][j][k] = 0.0;
			}
		}
	}

    //chain->NAA = 0; chain->aa = NULL; chain->erg = NULL; chain->xaa= NULL;
	//retv = getpdb(&(chain->aa), &(chain->NAA), &(chain->Nchains), infile);
    
    if(tempchain->NAA > 2){

	tempchain->aa[0].id = 'A';
	tempchain->aa[0].num = 0;
	tempchain->aa[0].chainid=0;
	tempchain->aa[0].etc=0;
	//fprintf(stderr,"PDBIN\n");
	//fprintf(stderr,"XAA\n");
	//for (int i = 1; i < tempchain->NAA; i++) {
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",tempchain->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");

	    freemem_chain(chain);
        allocmem_chain(chain,tempchain->NAA,tempchain->Nchains);

        copybetween(chain,tempchain);
        //fprintf(stderr,"%d amino acids in %d chains have been read in.\n", chain->NAA-1, chain->aa[chain->NAA-1].chainid);
        freemem_chain(tempchain);
        free(tempchain); 	

	    if (retv != EOF) {
// fixing peptide moved to main
//fprintf(stdout,"fixit = %d\n",sim_params->fixit);
//		    if(sim_params->fixit == 1) fixpeptide(chain->aa, chain->NAA, &(sim_params->protein_model));
//		    chkpeptide(chain->aa, chain->NAA, &(sim_params->protein_model));
		     return retv;
	    } else if (chain->NAA)		/* something was obtained earlier, move on */
		  return EOF;
	    else {
		    stop("Invalid input");
	    }
	}
   else{
//	fprintf(stderr, "pdbin: return EOF\n");
      freemem_chain(tempchain);
      free(tempchain); 
      return EOF;
    }

   /* we will never get here, but to avoid compiler warning... */
   return 0;
}

// mapping peptide onto the CRANKITE model: rebuilding the whole protein
// TODO: multi-chain proteins
void initialize(Chain *chain, Chaint * chaint, simulation_params *sim_params)
{
	int i;
	double rmse;
	vector orig;

//	fprintf(stderr,"INITIALIZE!\n");

	allocmem_chain(chain,chain->NAA,chain->Nchains);
    
	for (i = 1; i < chain->NAA; i++)
		fputc((chain->aa[i].id & COD) | (chain->aa[i].etc & ~COD), stderr);
	fprintf(stderr,"\n");

	/* before touching the protein, possibly save the sidechain dihedral angles,
	   because those cannot be calculated from only the backbone coordinates,
	   which we have when doing a crankshaft move */
	if ((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) {
	    for (i = 1; i < chain->NAA; i++) {
		initialise_sidechain_dihedral_angles((chain->aa)+i, &(sim_params->protein_model));
	    }
	}
	
	/* adjust Ca(i) atoms along the Ca(i-1)-Ca(i) and Ca(i)-Ca(i+1) directions
	   to 3.819 A.  Once that's done, also adjust the C and N distances, proportionally
	   to the ratio of the length of the Ca(i)-C(i)'s projection onto the Ca(i)-Ca(i+1) vector is. */
	/* this hardly changes the protein and any of its angles */
	rmse = repair_multichain(chain, chain);
	//fprintf(stderr,"CA repair RMSE %g",rmse);

	/* assign peptide bond orientations and adjust alpha-carbon locations */
	/* actually, this will rebuild the whole protein, keeping only the repaired Ca locations */
	amidorient(chain->xaa[0], NULL, (chain->aa) + 1);
	amidorient(chain->xaa_prev[1], NULL, (chain->aa) + 1);
	for (i = 1; i < chain->NAA - 1 ; i++) {
		chain->aa[i].SCRot = 0;
		if (chain->aa[i].chainid == chain->aa[i+1].chainid) { //build the next amino acid
			castvec(orig, chain->aa[i + 1].ca);
			//first find the right xaa[i]
			amidorient(chain->xaa[i], (chain->aa) + i, (chain->aa) + i + 1);
			//fprintf(stderr,"last chain.xaa %d %g %g %g\n",i,chain->xaa[i][0],chain->xaa[i][1],chain->xaa[i][2]);
			//then rebuild CA of aa[i+1]
			carbonate_f((chain->aa) + i + 1, (chain->aa) + i, chain->xaa[i]);
			//fprintf(stderr,"carb_f chain.ca %d %g %g %g\n",i+1,chain->aa[i+1].ca[0],chain->aa[i+1].ca[1],chain->aa[i+1].ca[2]);

			rmse += distance(orig, chain->aa[i + 1].ca);
			//fprintf(stderr,"rmse %g\n",rmse);
		} else { //for end of chain
			fprintf(stderr,"Orienting chain end differently! %d\n", i);
			amidorient(chain->xaa[i], (chain->aa) + i, NULL);
			fprintf(stderr,"Orienting chain start differently! %d\n", i+1);
			//leave the next chain's first CA (aa[i+1]) alone
			//and set up its xaa_prev orientation vector
			amidorient(chain->xaa_prev[chain->aa[i+1].chainid], NULL, (chain->aa) + i + 1);
		}
	}
	//i = chain->NAA - 1
	amidorient(chain->xaa[i], (chain->aa) + i, NULL);

	//fprintf(stderr,"XAA\n");
	//for (i = 1; i < chain->NAA; i++) {
	//	fprintf(stderr,"%d ",i);
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");

	/* input diagnostics: length, rmse */
	fprintf(stderr, " %d %g\n", chain->NAA - 1, sqrt(rmse / (chain->NAA - 1)));
	//fprintf(stderr, "CA RMSD of init %g\n", sqrt(rmse / (chain->NAA - 1)));

	/* build side chain */
	fulfill(chain,chaint, sim_params);
    
}

/* Fixed amino acids in list file will be marked with FIXED flag in aa.etc.
   These amino acids will not be moved in the MC simulation. */
void mark_fixed_aa_from_file(Chain *chain, simulation_params *sim_params) {

  //set all amino acids as not fixed (by default, the 0th amino acid is set as fixed.  TODO: check where it is set).
  for (int i=0; i<chain->NAA; i++) chain->aa[i].etc &= ~FIXED;

  if(!(sim_params->protein_model).fixed_aalist_file) return;

  fprintf(stderr,"marking fixed amino acids from file %s\n",(sim_params->protein_model).fixed_aalist_file);

  FILE *fptr = fopen((sim_params->protein_model).fixed_aalist_file, "r");
  if (!fptr) stop("mark_fixed_aa_from_file: problems while opening constraint file");

  fprintf(stderr, "Fixing amino acids:");
  int next;
  while (fscanf(fptr,"%d",&next) > 0) {
    if (next >= chain->NAA) stop("fixed amino acid beyond chain length");
    //fprintf(stderr,"constraining amino acid %d:",next);
    fprintf(stderr," %d%c",next,chain->aa[next].id);
    //fprintf(stderr,"  aa.etc before: %x",chain->aa[next].etc);
    chain->aa[next].etc |= FIXED;
    //fprintf(stderr,"  and after: %x\n",chain->aa[next].etc);
  }
  fprintf(stderr,"\n");
  fclose(fptr);

}


/* constrained amino acids in list file will be marked with CONSTRAINED flag in aa.etc */
void mark_constrained_aa_from_file(Chain *chain, simulation_params *sim_params) {

  if(!(sim_params->protein_model).external_constrained_aalist_file && !(sim_params->protein_model).external_constrained_aalist_file2) return;

  if ((sim_params->protein_model).external_constrained_aalist_file) {
    fprintf(stderr,"marking constrained amino acids from file %s\n",(sim_params->protein_model).external_constrained_aalist_file);

    FILE *fptr = fopen((sim_params->protein_model).external_constrained_aalist_file, "r");
    if (!fptr) stop("problems while opening constraint file");
  
    fprintf(stderr, "Constraining amino acids:");
    int next;
    while (fscanf(fptr,"%d",&next) > 0) {
      if (next >= chain->NAA) stop("constrained amino acid beyond chain length");
      //fprintf(stderr,"constraining amino acid %d:",next);
      fprintf(stderr," %d",next);
      //fprintf(stderr,"  aa.etc before: %x",chain->aa[next].etc);
      chain->aa[next].etc |= CONSTRAINED;
      //fprintf(stderr,"  and after: %x\n",chain->aa[next].etc);
    }
    fprintf(stderr,"\n");
    fclose(fptr);
  }
  if ((sim_params->protein_model).external_constrained_aalist_file2) {
    fprintf(stderr,"marking constrained amino acids from file %s\n",(sim_params->protein_model).external_constrained_aalist_file2);

    FILE *fptr = fopen((sim_params->protein_model).external_constrained_aalist_file2, "r");
    if (!fptr) stop("problems while opening constraint file");
  
    fprintf(stderr, "Constraining amino acids:");
    int next;
    while (fscanf(fptr,"%d",&next) > 0) {
      if (next >= chain->NAA) stop("constrained amino acid beyond chain length");
      //fprintf(stderr,"constraining amino acid %d:",next);
      fprintf(stderr," %d",next);
      //fprintf(stderr,"  aa.etc before: %x",chain->aa[next].etc);
      chain->aa[next].etc |= CONSTRAINED2;
      //fprintf(stderr,"  and after: %x\n",chain->aa[next].etc);
    }
    fprintf(stderr,"\n");
    fclose(fptr);
  }


}


void copybetween(Chain *to, Chain *from){
int i,j;
  if(to->NAA != from->NAA) {
	fprintf(stderr,"to->NAA=%d, from->NAA=%d\n",to->NAA,from->NAA);
	stop("Error in copybetween in metropolis.c\n");
  }
  if(to->Nchains != from->Nchains) {
	fprintf(stderr,"to->Nchains=%d, from->Nchains=%d\n",to->NAA,from->NAA);
	stop("Error in copybetween in metropolis.c\n");
  }
  for(j = 0; j < to->NAA; j++){
    to->aa[j].id = from->aa[j].id;
    to->aa[j].num = from->aa[j].num;
    to->aa[j].chainid = from->aa[j].chainid;
    to->aa[j].etc = from->aa[j].etc;
    to->aa[j].chi1 = from->aa[j].chi1;
    to->aa[j].chi2 = from->aa[j].chi2;
    for(i = 0; i < 3; i++){
      to->aa[j].h[i] = from->aa[j].h[i];	
	  to->aa[j].n[i] =  from->aa[j].n[i];		
	  to->aa[j].ca[i] = from->aa[j].ca[i];		
	  to->aa[j].c[i] = from->aa[j].c[i];		
	  to->aa[j].o[i] = from->aa[j].o[i];
	  to->aa[j].cb[i] = from->aa[j].cb[i];
	  to->aa[j].g[i] = from->aa[j].g[i];
	  to->aa[j].g2[i] = from->aa[j].g2[i];
	  to->xaa[j][i][0] = from->xaa[j][i][0];
	  to->xaa[j][i][1] = from->xaa[j][i][1];
	  to->xaa[j][i][2] = from->xaa[j][i][2];
    }		
  }		
  for(j = 0; j <= to->Nchains; j++){
      for(i = 0; i < 3; i++){
        to->xaa_prev[j][i][0] = from->xaa_prev[j][i][0];
        to->xaa_prev[j][i][1] = from->xaa_prev[j][i][1];
        to->xaa_prev[j][i][2] = from->xaa_prev[j][i][2];
      }
  }
  for(i = 0; i < to->NAA*to->NAA; i++){
	to->erg[i] = from->erg[i];  
  }
  to->ll = from->ll;	
}

