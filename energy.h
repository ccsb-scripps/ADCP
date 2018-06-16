/*
** This is an implementation of model interactions between two amino acids
** as well within a single amino acid. This is a rather simple force-field.
**
** Copyright (c) 2004 - 2008 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

/* CA-CA distance cutoff for vdW interactions */
//extern const double vdw_cutoff2_without_use_gamma;
//extern const double vdw_cutoff2_gg;
//extern const double vdw_cutoff2_gng;
/* CA-CA distance cutoff for Hbond interactions */
extern const double hbond_cutoff;


/* energy matrix and biasmap operations */
void energy_matrix_calculate(Chain *,Biasmap *, model_params *mod_params);
double totenergy(Chain *chain);
double locenergy(Chain *chain);
double extenergy(Chain *chain);
double targetenergy(Chain *chain);
double firstlastenergy(Chain *chain);

int bestRot(Chain *chain);
void energy_matrix_print(Chain *,Biasmap *, model_params *mod_params);
void biasmap_initialise(Chain *,Biasmap *, model_params *mod_params);
void biasmap_finalise(Biasmap *biasmap);
double centerX, centerY, centerZ, spacing;
//double *Nmapvalue, *Omapvalue, *CAmapvalue, *Hmapvalue, *Cmapvalue, *NAmapvalue, *Smapvalue, *emapvalue, *dmapvalue;
int NX, NY, NZ;
double targetBest, currTargetEnergy;
//double totalEBest;
double lowerGridEnergy(double);
void gridbox_initialise();
double *gridmapvalues[9];
void gridmap_initialise(char *, int);

double gridenergy(double X, double Y, double Z, int i, double charge);

void vectorProduct(float *a, float *b, float *c);
void normalizedVector(float *a, float *b, float *v);

float scoreSideChain(int nbRot, int nbAtoms, double *acharges, int *aTypes,  double coords[nbRot][nbAtoms][3], AA *a);

int getindex(int x, int y, int z);

/* wrappers for energy contributions at the amino acid level: energy contributions summed up */
/* energy of deformations and interactions within one amino acid */
double energy1(AA *, model_params *mod_params);
/* the energy of interactions between between two amino acids */
double energy2( Biasmap *,AA *,  AA *, model_params *mod_params);
/* the energy terms from terms that don't involve 1 or 2 residues */
double cyclic_energy(AA *, AA *, int);
double global_energy(int, int, Chain*, Chaint*,Biasmap *, model_params *mod_params);
double all_vdw(Biasmap *biasmap, Chain *chain, model_params *mod_params);

/* energy functions */
double linear_decay(double, double, double); 
double stress(AA *a, model_params *mod_params);
double proline(AA *a, AA *b);
double bias(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params);
double hbond(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params);
double hydrophobic(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params);
int hydrophobic_interaction_intensity(AA *a, AA *b, model_params *mod_params);
double hydrophobic_low(double distance, double contact_cutoff, model_params *mod_params);
double electrostatic(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params);
double sidechain_hbond(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params);
double secondary_radius_of_gyration(int start, int end, Chain *chain, Chaint *chaint, Biasmap *biasmap, model_params *mod_params, int which_atom, int only_hydrophobic);

/* custom tests that depend on the energy implementation */
void energy_probe_1(Chain *chain,Biasmap *biasmap,simulation_params *sim_params);
void energy_contributions_in_energy_c(Chain * chain,Biasmap *biasmap, double tote, model_params *mod_params, FILE *outfile);
void exclude_energy_contributions_in_energy_c(Chain * chain,Biasmap *biasmap, double tote, model_params *mod_params, FILE *outfile);
