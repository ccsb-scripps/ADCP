/*
** ML inference of model parameters using Contrastive Divergence.
**
** Copyright (c) 2004 - 2008 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

typedef struct {

  int iter_max;
  char *learn_string;
  char *restart_filename;
  int iter_start;
  //int restart_frequency;
  char *pdb_list_filename;
  char *output_filename;

  double aa[36]; //param_value[36]; //aa
  double raa[36]; //param_delta[36]; //ra
  double total_energy_gradient[36];
  double previous_move[36];
  //int n_proteins;
  double adaptive_learning_param;

} cdlearn_params;

void update_sim_params_from_cd_learn(simulation_params *sim_params, cdlearn_params *cd_params);
void cd_learn_param_initialise(cdlearn_params *this);
void cd_learn_write_restart_file(cdlearn_params *this, char *restart_filename);
void cd_param_finalise(cdlearn_params *this);
void print_contact_map(double *distb, int NAA, FILE *outfile);
