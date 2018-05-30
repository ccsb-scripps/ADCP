/*
** This contains the description of polypeptide structure
** and IO routines.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

/* constants to store properties of the amino acid */
#define COD 0x1F /* amino acid code */
#define PSI 0x20 /* if the amino acid is in not helix conformation (psi dihedral angle is not in [-120deg,60deg]) */
#define LEV 0x40 /* L amino acid */
#define CIS 0x80 /* cis peptide bond */

#define HYDROPHOBIC   0x10000 /* hydrophobic side chain */
#define AMPHIPATHIC   0x20000 /* amphipathic side chain */
#define ELECTROSTATIC 0x40000 /* charged side chain     */
#define CONSTRAINED   0x80000 /* constrained amino acid (see external potential)   */
#define CONSTRAINED2  0x100000 /* constrained amino acid (see external potential)   */
#define FIXED         0x200000 /* fixed amino acid */

/* fixed Ca-Cb bond length and Ca-Ca distance */
extern const double cacb, caca, caca_;
/* fixed trans peptide bond atom coordinates in local frames */
extern vector ca, h, n, c, o;
/* fixed cis (proline) peptide bond atom coordinates */
extern vector ca_, h_, n_;
/* normalized bond vectors */
extern vector n1, c1, n1_, cn1, cn1_;
/* approximate reciprocal Ca-N and Ca-C bond lengths and skew between them */
extern double ican, icac, ican_, skew, skew_;


/* amino acid type */
typedef struct _AA {
	vector h;		/* hydrogen */
	vector n;		/* nitrogen */
	vector ca;		/* alpha carbon */
	vector c;		/* carbonyl carbon */
	vector o;		/* carbonyl oxygen */
	vector cb;		/* beta carbon */
	vector g;		/* gamma carbon, oxygen, or sulfur */
	vector g2;		/* second gamma atom
				   In case there are 2 (for ILE, THR and VAL)
				   ILE: g: ethyl group (CG1), g2: methyl group (CG2)
				   THR: g: oxygen (OG1), g2: carbon (CG2)
				   VAL: symmetric (CG1 and CG2) */
	double chi1;		/* n-ca-cb-g side chain dihedral angle */
	double chi2;		/* n-ca-cb-g2 side chain dihedral angle */
	int etc;		/* cis-trans, levo-dextro, etc flags */
	int num;		/* sequence position */
	char id;		/* 1-letter type abbreviation */
	int chainid;		/* chain id (1, 2, ...) */
	int SCRot;   /*best scoring side chain coords*/
	int donor; /*H bond donor, -N-H-----O*/
	int acceptor; /*H bond accptor -C=O-----H*/
} AA;

typedef struct _FLEX_data{
  int* oxy_index;
  int *Hbond_aaH, *Hbond_aaO;
  int number_hbond;
  int total_flex;
  int accepted_flex;
  int read_in_flex;
} FLEX_data;


/* amino acid chain type */
typedef struct _Chain {
	AA *aa;
	triplet *xaa;
	triplet *xaa_prev; //previous xaa for chain start only
	double *erg;
	double ll; /*logL only used for Nested sampling */
    int NAA;
    int Nchains;
    FLEX_data *flex_data; //only used in nma.c otherwise ignored
} Chain;

/* temporary amino acid chain type */
typedef struct _Chaint {
	triplet *xaat;
	triplet *xaat_prev; //previous xaa for chain start only
	AA *aat;
	double *ergt;
} Chaint;

/* bias map type */
typedef struct _Biasmap {
	double *distb;
    int NAA;
} Biasmap;


/* initialisation of constants */
void peptide_init();

/* geometry tests */
void chkssbond( AA *a, int count);
void chkpeptide( AA *, int, model_params *mod_params);

/* amino acid and peptide building and modification */
void carbonate_f( AA *,  AA *, triplet);
void carbonate_b( AA *,  AA *, triplet);
void acidate( AA *, triplet, triplet, simulation_params *sim_params);
void gammalate( AA *a, double chi, int which_gamma, model_params *mod_params);
void fixpeptide( AA *, int, model_params *mod_params);
void adjust_to_closest_state(double *chi, char aa_type);
void initialise_sidechain_dihedral_angles( AA *a, model_params *mod_params);
void ramaset(triplet, double, double, double, int);
void omegaset(triplet, double, int);
void amidorient(triplet,  AA *,  AA *);
void aat_init(Chain *chain,Chaint *chaint);
void fulfill(Chain *chain, Chaint *chaint, simulation_params *sim_params);
void build_peptide_from_sequence(Chain *chain, Chaint *chaint, char *str, simulation_params *sim_params);

/* setup and i/o */
void initialize(Chain *chain, Chaint *chaint, simulation_params *sim_params);
int pdbin(Chain *chain, simulation_params *sim_params, FILE *infile);
void allocmem_chain( Chain *chain, int NAA, int Nchains);
void freemem_chain(Chain *chain);
void freemem_chaint(Chaint *chaint);
void mark_fixed_aa_from_file(Chain *chain, simulation_params *sim_params);
void mark_constrained_aa_from_file(Chain *chain, simulation_params *sim_params);
void copybetween(Chain *,Chain*);

/* more geometry tests */
int hdonor( AA *,  AA *, model_params *mod_params);
double hstrength(vector,vector,vector,vector, model_params *mod_params);
int contact( AA *,  AA *, model_params *mod_params);
int aligned( AA *,  AA *);

/* peptide chain i/o */
int pdbrecord( AA *, int, model_params *mod_params, FILE *outfile);
void pdbprint( AA *, int, model_params *mod_params, FILE *outfile, double *totenergy);
int getaa( AA *, FILE *infile);
int getpdb( AA **, int *NAA, int *Nchains, FILE *infile);
