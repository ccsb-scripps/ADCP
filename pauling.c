/*
**  This program builds a PDB model from Ramachandran and other angles.
**  When ramachandran angles are omitted, a random conformation is produced.
**  The input lines should begin from one-letter amino acid code followed by
**  phi, psi, tau, omega, and chi angles in degrees.
**
**  Lipa 1.2, Copyright (c) 2004-2010 Alexei Podtelezhnikov
*/

#define VER "Lipa 1.2, Copyright (c) 2004-2010 Alexei Podtelezhnikov\n"
#define USE "Usage: %s [[-f] filein] [-o fileout] [-s seed] [-m model]\n"

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<ctype.h>

#include"params.h"
#include"aadict.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"

unsigned int seed = 0, model = 1;
double pool[4] = { 2 * M_PI_3, -2 * M_PI_3, 2 * M_PI_3 };
double NaN;

void read_options(int argc, char *argv[], int *use_gamma_atoms, int *use_rand_tau)
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
		case 'm':
			sscanf(argv[i], "%u", &model);
			break;
		case 's':
			sscanf(argv[i], "%u", &seed);
			break;
		case 'g':
			//sscanf(argv[i],"%d",use_gamma_atoms);
			if (strcmp(argv[i],"NONE")==0) {
				*use_gamma_atoms = NO_GAMMA;
			} else if (strcmp(argv[i],"LINUS_GAMMA")==0) {
				*use_gamma_atoms = LINUS_GAMMA;
			} else if (strcmp(argv[i],"CORRECT_GAMMA")==0) {
				*use_gamma_atoms = CORRECT_GAMMA;
			} else if (strcmp(argv[i],"CORRECT_KMQR_GAMMA")==0) {
				*use_gamma_atoms = CORRECT_KMQR_GAMMA;
			}
			if (*use_gamma_atoms==NO_GAMMA) {
				fprintf(stderr,"WARNING! Not using gamma atoms.\n");
			}
			break;
		case 't':
			sscanf(argv[i],"%d",use_rand_tau);
			break;
		default:
			fprintf(stderr, VER USE, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
}

static double randphi(void)
{
	/*For nested sampling, need genuinely random dihedral angles */
	return (((double)rand() / RAND_MAX)-0.5)*360;
	/*return 120 * M_PI_180 * rand() / RAND_MAX - 180 * M_PI_180; */
}

static double randtau(void)
{
	/*For nested sampling, need genuinely random valence angles */
	return (((double)rand() / RAND_MAX)-0.5)*360;
	/* for energy minimum use */ 
	//return atan2(M_SQRT2, -0.5);
}

static double randpsi(void)
{
	/*For nested sampling, need genuinely random dihedral angles */
	return (((double)rand()/RAND_MAX) - 0.5)*360;
	/*double rpsi, s4;

	rpsi = 200 * M_PI_180 * rand() / RAND_MAX;

	if (rpsi < 100 * M_PI_180) {
		s4 = pool[0] + pool[1] + pool[2] + pool[3];
		if (-410 * M_PI_180 < s4 && s4 < -280 * M_PI_180)
			rpsi += 100 * M_PI_180;
		else
			rpsi -= 80 * M_PI_180;
	}

	return rpsi; */
}

static double randchi(char id, sidechain_properties_ *sidechain_properties)
{ return sidechain_dihedral(id, sidechain_properties);
}

void pdbrecur(char type, double phi, double psi,
	      double tau, double omega, double chi, double chi2, simulation_params *sim_params)
{
	matrix t;
	triplet x0, xn;
	static int j = 0;
	static triplet xc = { {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.} };
	static AA a = {.ca = {0., 0., 0.},.num = 0,.etc = LEV };

	if (isnan(omega))
		omega = M_PI;
	if (isnan(phi))
		phi = 1.0;
	if (isnan(psi))
		psi = -4.0;

	a.id = type;
	a.num++;

	if (-M_PI / 2 < omega && omega < M_PI / 2)
		a.etc |= CIS;
	else
		a.etc &= ~CIS;

	omegaset(x0, omega, a.etc & CIS);
	transset(t, xc);
	rotation(xn, t, x0);

	/* a little bit of maigic that implements the omega dihedral */
	lincomb(x0[1], 0.5223, xc[1], 0.4777, xn[1]);
	lincomb(x0[2], 0.5223, xc[2], 0.4777, xn[2]);
	carbonate_f(&a, &a, x0);

	ramaset(x0, phi, psi, tau, a.etc & CIS);
	transset(t, xn);
	rotation(xc, t, x0);

	/* set backbone and sidechain atoms (incl. gamma, default) */
	acidate(&a, xn, xc, sim_params);
	a.etc |= CA_ | CB_ | N__ | O__ | C__;
	/* update gamma atoms */
	if (isnan(chi))
		a.g[2] = NaN;	/* invalidate */
	else {
		gammalate(&a, chi, 1, &(sim_params->protein_model));
		a.etc |= G__;
		if (isnan(chi2))
			a.g2[2] = NaN;	/* invalidate */
		else {
			gammalate(&a, chi2, 2, &(sim_params->protein_model));
			a.etc |= G2_;
		}
	}

	j = pdbrecord(&a, j, &(sim_params->protein_model), stdout);
	if (psi == -4.0)
		printf("TER   %5d      %3s A%4d\n", ++j, aa123(a.id), a.num);

	/* fix accumulating errors every 32 residues */
	if ((a.num & 0x1F) == 0)
		fixtriplet(xc);
	/* printout(xc); */
}

char uppercase(char this) {
  if (isupper(this)) return this;
  if (islower(this)) return toupper(this);
  fprintf(stderr,"Got non-alpha character >%c<.",this);
  exit(EXIT_FAILURE);

}

int main(int argc, char *argv[])
{
	
	
	int count;
	char line[1024], id;
	double phi, psi, tau, omega, chi, chi2, tau0;
	const double coef = M_PI / 180.0;
	int use_rand_tau;
	simulation_params sim_params;
    int use_gamma_atoms = LINUS_GAMMA;
	use_rand_tau = 0;
	
	read_options(argc, argv, &use_gamma_atoms, &use_rand_tau);

	param_initialise(&sim_params);
	sim_params.protein_model.use_gamma_atoms = use_gamma_atoms;
	initialize_sidechain_properties(&(sim_params.protein_model));
	
	if (seed == 0)
		seed = (unsigned int) time(NULL);
	srand(seed);
    
	peptide_init();

	tau0 = atan2(M_SQRT2, -0.5);
	NaN = strtod("NaN", NULL);

	printf("MODEL     %4d\n", model);
	while (fgets(line, sizeof(line), stdin) != NULL) {
		count = sscanf(line, "%c %lg %lg %lg %lg %lg %lg",
			       &id, &phi, &psi, &tau, &omega, &chi, &chi2);

		switch (count) {	/* convert to radians */
		case 7:
			chi2 *= coef;
		case 6:
			chi *= coef;
		case 5:
			omega *= coef;
		case 4:
			tau *= coef;
		case 3:
			psi *= coef;
		case 2:
			phi *= coef;
			break;
		case 1:
			if (line[1] == '\n' || line[1] == '\r')
				break;
			if (line[0] == 'I' && line[1] == 'd')
				continue; /* skip header */
			fprintf(stderr, "Invalid input\n");
			return EXIT_FAILURE;
		default:
			continue;	/* read next line */
		}

		switch (count) {	/* fill the blanks */
		case 1:
			phi = randphi();
		case 2:
			pool[3] = phi;
			psi = randpsi();
		case 3:
			if(use_rand_tau) tau = randtau();
			else tau = tau0;
		case 4:
			if (islower(id)) {
			  omega = 0;
			} else {
			  omega = M_PI;
			}
		case 5:
			if(use_gamma_atoms != NO_GAMMA){
			  chi = randchi(id,sim_params.protein_model.sidechain_properties);
			}
			else{
			  chi = NaN;
		    }
		case 6:
			if(use_gamma_atoms != NO_GAMMA){
			  if(id == 'V' || id == 'T' || id == 'I'){
			    chi2 = sidechain_dihedral2(id,chi, sim_params.protein_model.sidechain_properties);
			  }
			  else chi2 = NaN;  
			}
			else{
			  chi2 = NaN;
		    }
		}

		pdbrecur(uppercase(id), phi, psi, tau, omega, chi, chi2, &sim_params);

		pool[0] = pool[2];
		pool[1] = phi;
		pool[2] = psi;
	}
	printf("ENDMDL\n");
	return EXIT_SUCCESS;
}
