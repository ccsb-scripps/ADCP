/*
** This is an implementation of model interactions between two amino acids
** as well within a single amino acid. This is a rather simple force-field.
**
** Copyright (c) 2004 - 2008 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

extern const int skip_14_vdw;
/* energy penalty for clashing atoms */
extern const double clash_energy_at_hard_cutoff;


/* wrappers for energy contributions at the amino acid level: energy contributions summed up */
//double all_vdw(Biasmap *biasmap, Chain *chain, model_params *mod_params);

/* energy functions */
double clash(AA *a, model_params *mod_params);
double vdw_hard_cutoff(vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff);
double vdw_lj(vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff);
double exclude_neighbor(AA *a, AA *b, model_params *mod_params);
double exclude(AA *a, AA *b, double d2, model_params *mod_params);

/* custom tests that depend on the energy implementation */
//void exclude_energy_contributions_in_energy_c(Chain * chain,Biasmap *biasmap, double tote, model_params *mod_params, FILE *outfile);

/* vdw cutoff tests */
double vdw_low_level(vector CA_1, vector CA_2, vector A_1, vector A_2, double radii_1, double radii_2);
double vdw_backbone_constants(Chain* chain, model_params *mod_params, FILE *outfile, int verbose);
double vdw_gamma_gamma(AA *a, AA *b, model_params *mod_params);
double vdw_gamma_nongamma( AA *a,  AA *b, model_params *mod_params);
void vdw_maxgamma_calc(Chain *chain, model_params *mod_params, FILE *outfile, int verbose);
void vdw_cutoff_distances_calculate(simulation_params *sim_params, FILE *outfile, int verbose);
void update_sim_params_from_chain(Chain *chain,simulation_params *sim_params);
#ifdef LJ_HBONDED_HARD
double exclude_hard(AA *a, AA *b, double d2, model_params *mod_params, int hbond_proximity);
#endif
