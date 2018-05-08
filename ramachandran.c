/*
**  This program calculates Ramachandran angles for a given PDB file.
**  Parsing options and default values are as follows.
**
**  Rama 1.2, Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
*/

#define VER "Rama 1.2, Copyright (c) 2004 - 2010 Alexei Podtelezhnikov\n"
#define USE "Usage: %s [[-f] filein] [-o fileout]\n"

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"params.h"
#include"aadict.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"

double tetrahedron(vector a0, vector a1, vector a2)
{
	vector b0, b1;

	subtract(b0, a0, a1);
	subtract(b1, a2, a1);

	return M_180_PI * angle(b0, b1);
}

double ramachandran(vector a0, vector a1, vector a2, vector a3)
{
	return M_180_PI * dihedral_4(a0, a1, a2, a3);
}

double solid(vector a0, vector a1, vector a2, vector a3)
{
	vector b1, b2, b3;

	subtract(b1, a1, a0);
	subtract(b2, a2, a0);
	subtract(b3, a3, a0);

	return M_1_PI * excess(b1, b2, b3);
}

double *normal(vector n, vector a0, vector a1, vector a2)
{
	vector b1, b2;

	subtract(b1, a1, a0);
	subtract(b2, a2, a0);

	crossprod(n, b1, b2);

	return n;
}

int iscis(vector a0, vector a1, vector a2, vector a3)
{
	vector n1, n2;

	normal(n1, a0, a1, a2);
	normal(n2, a1, a2, a3);

	return dotprod(n1, n2) > 0;
}

/* ramafix can be used to analyze accumulating errors in orientations */
void ramafix(AA *a,  AA *b,
	     double *phi_, double *psi_, double *tau_, double *omega_)
{
	matrix t;
	triplet x0, xn, xr;
	static int i = 0;
	static triplet xc;
	double phi, psi, tau, omega;

	if (iscis(a->ca, a->c, b->n, b->ca))
		b->etc |= CIS;
	else
		b->etc &= ~CIS;
	amidorient(xr, a, b);
	/* printout(xr); */

	if (i++ == 0) {
		casttriplet(xc, xr);
		peptide_init();
		return;
	}

	phi = *phi_ * M_PI_180;
	psi = *psi_ * M_PI_180;
	tau = *tau_ * M_PI_180;
	omega = *omega_ * M_PI_180;

	omegaset(x0, omega, a->etc & CIS);
	transset(t, xc);
	rotation(xn, t, x0);

	ramaset(x0, phi, psi, tau, a->etc & CIS);
	transset(t, xn);
	rotation(xc, t, x0);

	/* fix accumulating errors every 32 residues */
	if ((i & 0x1F) == 0)
		fixtriplet(xc);
	/* printout(xc); */

	fprintf(stderr, "Twist = %g    Bend = %g\n",
		euler_twist(xr, xc), euler_bend(xr, xc));
}

void read_options(int argc, char *argv[])
{
	int i, opt;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			freopen(argv[i], "r", stdin);
			continue;
		}

		opt = argv[i][1];
		if (++i >= argc)
			opt = 0;

		switch (opt) {
		case 'f':
			freopen(argv[i], "r", stdin);
			break;
		case 'o':
			freopen(argv[i], "w", stdout);
			break;
		default:
			fprintf(stderr, VER USE, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char *argv[])
{
	int status_a, status_b;
	AA a, b;
	double phi, psi, tau, omega, chi,chi2;
	double caca, cbcb, sigma, eta, baab, naac, sol;
	vector n1, n2;
	double NaN;
	const char *frmt =
	    "%c%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f%7.3f%7.1f%7.3f%7.3f%7.1f%7.1f\n";

	read_options(argc, argv);

	/* NaN means that an atom is missing */
	NaN = strtod("NaN", NULL);

	do
		status_b = getaa(&b,stdin);
	while (!(b.etc & CA_) && status_b != EOF);

	if (b.etc & CA_) {
		printf("Id   phi    psi    tau  omega    chi    chi2  solid    eta"
		       "  CA-CA  CB-CB  BA^AB  NA^AC\n");
		goto start;
	} else {
		fprintf(stderr, "No amino acids recognized\n");
		return EXIT_SUCCESS;
	}

	while (1) {
		tau = tetrahedron(b.n, b.ca, b.c);
		sol = solid(b.ca, b.n, b.c, b.cb); //Had to remove distance here
 		if (b.id == 'G' || b.id == 'A')  //  || distance(b.cb, b.g) > 4.0))
			chi = NaN;
		else		/* gauche+ corresponds to +60 degrees (IUPAC) */
			chi = ramachandran(b.n, b.ca, b.cb, b.g);
		if(b.id == 'T' || b.id == 'V' || b.id == 'I')
			chi2 = ramachandran(b.n, b.ca, b.cb, b.g2);
		else 
		    chi2 = NaN;

		a = b;
		status_a = status_b;
		do
			status_b = getaa(&b,stdin);
		while (!(b.etc & CA_) && status_b != EOF);
		if (!(b.etc & CA_) && status_b == EOF)
			status_a = EOF;

		if (status_a != EOF && status_a != 0x7FFF && status_a != 0x7FFE) {
			psi = ramachandran(a.n, a.ca, a.c, b.n);

			normal(n2, a.ca, a.c, b.n);
			eta = M_180_PI * angle(n1, n2);

			/* ramafix(&a, &b, &phi, &psi, &tau, &omega); */

			printf(frmt, a.id, phi, psi, tau, omega, chi, chi2, sol, eta,
			       caca, cbcb, baab, naac);
			
			normal(n1, b.ca, b.n, a.c);
			omega = ramachandran(a.ca, a.c, b.n, b.ca);
			phi = ramachandran(a.c, b.n, b.ca, b.c);
			//fprintf(stderr,"(1) %g\n",phi);
			sigma = psi + omega + phi;

			if (sigma > 180.0)
				sigma -= 360.0;
			if (sigma < -180.0)
				sigma += 360.0;

			if (b.id == 'G' || distance(b.ca, b.cb) > 4.0)
				b.cb[0] = NaN;

			baab = ramachandran(a.cb, a.ca, b.ca, b.cb);
			naac = ramachandran(a.n, a.ca, b.ca, b.c);
			caca = sqrt(distance(a.ca, b.ca));
			cbcb = sqrt(distance(a.cb, b.cb));
		} else {
			psi = eta = NaN;

			printf(frmt, a.id, phi, psi, tau, omega, chi, chi2, sol, eta,
			       caca, cbcb, baab, naac);

			if (status_a == EOF)
				break;

		      start:
			omega = phi = sigma = n1[0] = chi = chi2 = NaN;
			//fprintf(stderr,"(2) %g\n",phi);
			baab = naac = caca = cbcb = NaN;
		}
	}

	return EXIT_SUCCESS;
}
