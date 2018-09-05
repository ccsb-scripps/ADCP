/*
** Amino acid lexicon and conversions.
**
** Copyright (c) 2004 - 2007 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

/* 3-letter amino acid codes */
extern const char *aac[];
extern const char aa_id[]; 
extern const char *unu[];
#define aa123(x) aac[(x) & 0x1F]

/* atom binary constants */
#define ATM 0xFF00
#define H__ 0x100
#define N__ 0x200
#define C__ 0x400
#define O__ 0x800
#define CA_ 0x1000
#define CB_ 0x2000
#define G__ 0x4000
#define G2_ 0x8000

#define BACKBONE_ACCEPTOR_RADIUS 2.25
#define BACKBONE_DONOR_RADIUS 2.25


/* amino acid code conversion */
char aa321(char *cod);
int convert_to_index(char id);

/* sidechain property library initialisation */
void initialize_one_sidechain_properties(sidechain_properties_ *sidechain_properties,
					 int num, char id, double number_of_gamma,
					 double beta_gamma1_distance,  double alpha_beta_gamma1_angle,
					 double beta_gamma2_distance,  double alpha_beta_gamma2_angle,
					 double sidechain_vdw_radius1,  double sidechain_vdw_radius2,
					 double sidechain_vdw_depth1,  double sidechain_vdw_depth2,
					 double sidechain_dihedral_gauche_plus_prob, double sidechain_dihedral_gauche_minus_prob,
					 int hydrophobic_atoms, double hydrophobic_contact_radius_CB, double hydrophobic_contact_radius_G1, double hydrophobic_contact_radius_G2,
					 int charged_atom, double charged_atom_charge,
					 int hydrogen_bond_donor_atoms, int hydrogen_bond_acceptor_atoms,
					 double hydrogen_bond_donor_radius, double hydrogen_bond_acceptor_radius );
void initialize_sidechain_properties(model_params *mod_params);

/* sidechain property query functions */
double charge(char id, sidechain_properties_ *sidechain_properties);
int beta_gamma_dist(char id, int which_gamma, double *r, double *theta, sidechain_properties_ *sidechain_properties);
double sidechain_dihedral(char id, sidechain_properties_ *sidechain_properties);
double sidechain_dihedral2(char id, double chi, sidechain_properties_ *sidechain_properties);
double sidechain_vdw_radius(char id, int which_gamma, sidechain_properties_ *sidechain_properties);
double sidechain_vdw_depth(char id, int which_gamma, sidechain_properties_ *sidechain_properties);
double sidechain_vdw_depth_sqrt(char id, int which_gamma, sidechain_properties_ *sidechain_properties);
double hydrophobic_contact_radius(char id, int atom, sidechain_properties_ *sidechain_properties);
unsigned int hydrophobic_atoms_list(char id, sidechain_properties_ *sidechain_properties);
int hbond_donor(char id, int atom, sidechain_properties_ *sidechain_properties);
int hbond_acceptor(char id, int atom, sidechain_properties_ *sidechain_properties);
double sidechain_hbond_donor_radius(char id, sidechain_properties_ *sidechain_properties);
double sidechain_hbond_acceptor_radius(char id, sidechain_properties_ *sidechain_properties);
