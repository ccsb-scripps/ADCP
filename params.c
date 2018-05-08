/*
** Module with all the simulation parameters **
** and its default values
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include<float.h>
#include"error.h"
#include"params.h"



/***********************************************************/
/****       INITIALISING AND FINALISING  ROUTINES       ****/
/***********************************************************/

/* initialising simulation parameters (including model parameters) */
/* default values should be set here */
void param_initialise(simulation_params *this) {
  /* general simulation */
  this->infile = stdin;
  this->outfile = stdout;
  this->infile_name = NULL;
  this->outfile_name = NULL;
  this->pace = 0;
  this->stretch = 16;
  this->tmask = 0x0;
  this->seed = 0;
  this->prm = NULL;
  this->acceptance_rate = 0.5;
  this->amplitude = -0.1;
  this->keep_amplitude_fixed = 0;
  this->acceptance = 1.0;
  this->accept_counter = 0;
  this->reject_counter = 0;
  this->acceptance_rate_tolerance = 0.03;
  this->amplitude_changing_factor = 0.9;
  this->number_initial_MC = 0;
  this->MC_lookup_table = NULL; //lookup table for random moves
  this->MC_lookup_table_n = NULL; //number of valid elements in the lookup table for random moves
  /* peptide */
  this->seq = NULL;
  this->sequence = NULL;
  this->NAA = 0;
  this->Nchains = 0;
  /* thermodynamic beta */
  this->thermobeta = 0;
  this->lowtemp = 0;
  this->beta1 = 1;
  this->beta2 = 0;
  this->bstp = 1;
  this->intrvl = 16384;
  this->energy_gradient = malloc(36*sizeof(double));
  this->energy_probe_1_this = malloc(36*sizeof(double));
  this->energy_probe_1_last = malloc(36*sizeof(double));
  this->energy_probe_1_calc = malloc(36*sizeof(int));
  int i;
  for (i=0; i<36; i++) {
    this->energy_gradient[i] = 0;
    this->energy_probe_1_this[i] = 0;
    this->energy_probe_1_last[i] = 0;
    this->energy_probe_1_calc[i] = 1;
  }
  /* nested sampling */
  this->NS = 0;
  this->iter = 0;
  this->iter_start = 0;
  this->iter_max = 1000;
  this->logLstar = DBL_MAX;
  this->logZ = -DBL_MAX;
  this->logfactor = 1;
  this->alpha = 1;
  this->logX = 0;
  this->logX_start = 0;
  this->Delta_logX = 0;
  this->log_DeltaX = 0;
  this->H = 0;
  this->N = 0;

  /* checkpointing */
  this->num_NS_per_checkpoint = 0;
  this->checkpoint_filename = NULL;
  this->checkpoint_file = NULL;
  this->checkpoint_counter = 0;
  this->restart_from_checkpoint = 0;
  this->checkpoint = 0;

  model_param_initialise(&(this->protein_model));
  flex_params_initialise(&(this->flex_params));

}

/* finalising simulation parameters (including model parameters) */
/* default values should be set here */
void param_finalise(simulation_params *this) {
  /* general simulation */

  if (this->infile && this->infile!=stdin) {
    fclose(this->infile);
    this->infile = NULL;
  }

  if (this->outfile && this->outfile!=stdout) {
    fclose(this->outfile);
    this->outfile = NULL;
  }

  if (this->infile_name) free(this->infile_name);
  if (this->outfile_name) free(this->outfile_name);

  this->pace = 0;
  this->stretch = 0;
  this->tmask = 0x0;
  this->seed = 0;

  if (this->prm) free(this->prm);

  this->acceptance_rate = 0.;
  this->amplitude = 0.;
  this->keep_amplitude_fixed = 0;
  this->acceptance = 0.;
  this->accept_counter = 0;
  this->reject_counter = 0;
  this->acceptance_rate_tolerance = 0;
  this->amplitude_changing_factor = 0;
  this->number_initial_MC = 0;
  if (this->MC_lookup_table) free(this->MC_lookup_table);
  if (this->MC_lookup_table_n) free(this->MC_lookup_table_n);


  /* peptide */

  if (this->seq) free(this->seq);
  if (this->sequence) free(this->sequence);
  this->NAA = 0;
  this->Nchains = 0;
  /* thermodynamic beta */
  this->thermobeta = 0;
  this->lowtemp = 0;
  this->beta1 = 0;
  this->beta2 = 0;
  this->bstp = 0;
  this->intrvl = 0;
  if (this->energy_gradient) free(this->energy_gradient);
  if (this->energy_probe_1_calc) free(this->energy_probe_1_calc);
  if (this->energy_probe_1_this) free(this->energy_probe_1_this);
  if (this->energy_probe_1_last) free(this->energy_probe_1_last);
  /* nested sampling */
  this->NS = 0;
  this->iter = 0;
  this->iter_start = 0;
  this->iter_max = 0;
  this->logLstar = 0;
  this->logZ = 0;
  this->logfactor = 0;
  this->alpha = 0;
  this->logX = 0;
  this->logX_start = 0;
  this->Delta_logX = 0;
  this->log_DeltaX = 0;
  this->H = 0;
  this->N = 0;

  /* checkpointing */
  this->num_NS_per_checkpoint = 0;
  if (this->checkpoint_filename) free(this->checkpoint_filename);
  if (this->checkpoint_file) {
	fclose(this->checkpoint_file);
	this->checkpoint_file = NULL;
  }
  this->checkpoint_counter = 0;
  this->restart_from_checkpoint = 0;
  this->checkpoint = 0;


  model_param_finalise(&(this->protein_model));
  flex_params_finalise(&(this->flex_params));

}


/* initialise model parameters */
void model_param_initialise(model_params *this) {

  /* gamma atoms */
  this->fixit = 1;
  this->use_gamma_atoms = LINUS_GAMMA;
  this->use_original_gamma_atoms = 1;
  this->use_3_states = 1;
  this->fix_chi_angles = 0;
  this->fix_CA_atoms = 0;
  /* atomic radii */
  /* radii from Ward et al, 1999 */
  // rca = 1.75, rcb = 1.75, rc = 1.65, rn = 1.55, ro = 1.40; 
  /* radii from Hopfinger, 1973 */
  // rca = 1.57, rcb = 1.57, rc = 1.42, rn = 1.29, ro = 1.29
  // rs = 1.8 from wikipedia I'm sorry ,NB
  /* LINUS */
  // rca = 1.85, rcb = 2.0, rc = 1.85, rn = 1.75, ro = 1.6, rs = 2.0; 
  this->rca = 2.43;
  this->rcb = 1.97;
  this->rc = 1.82;
  this->rn = 1.74;
  this->ro = 1.98;
  this->rs = 3.10;
  this->rring = 2.00;
  this->vdw_depth_ca = 0.018;
  this->vdw_depth_cb = 0.018;
  this->vdw_depth_c = 0.018;
  this->vdw_depth_n = 0.018;
  this->vdw_depth_o = 0.018;
  this->vdw_depth_s = 0.018;
  this->vdw_depth_ring = 0.2;

  /* The vdW cutoff distances have to be calculated later, because
     the cutoff calculating routine depends on vdw.c */
  this->vdw_gamma_gamma_cutoff = NULL;
  this->vdw_gamma_nongamma_cutoff = NULL;
  this->vdw_backbone_cutoff = 500;
  this->vdw_use_extended_cutoff = 0;
  this->vdw_extended_cutoff = DEFAULT_EXTENDED_VDW_CUTOFF_GAMMA;
  this->vdw_potential = LJ_VDW_POTENTIAL;
  this->vdw_function = NULL;
  this->vdw_clash_energy_at_hard_cutoff = 30; //default value for LJ
  this->vdw_lj_neighbour_hard = 0;
  this->vdw_lj_hbonded_hard = 0;
  
  this->rel_vdw_cutoff = 2.0;
  this->vdw_uniform_depth = 0;
  vdw_param_zero(this);
  vdw_param_calculate(this);

  /* stress */
  this->stress_k = 98.;
  this->stress_angle = 1.20427718387608740808;

  /* hydrogen bond */
  this->hboh2 = 4.04;
  this->hbohn = 0.928;
  this->hbcoh = 0.772;
  //this->hboh_decay_width = 1.0;
  //this->hbohn_decay_width = 0.1;
  //this->hbcoh_decay_width = 0.1;
  this->hbs = 4.98;
  /* contact parameters */
  this->touch2 = 38.44;
  this->part = 1.0;
  this->split = -1.0;
  this->sts = 0.0;
  /* biasing force constants */
  this->contact_map_file = NULL;
  this->bias_eta_beta = 3.7;
  this->bias_eta_alpha = 15.3;
  this->bias_kappa_alpha_3 = 0.0;
  this->bias_kappa_alpha_4 = 0.0;
  this->bias_kappa_beta = 0.85;
  this->prt = 1.;
  this->bias_r_alpha = 5.39;
  this->bias_r_beta = 5.39;
  /* hydrophobicity */
  this->kauzmann_param = 0.122;
  this->hydrophobic_cutoff_range = 2.8;
  this->hydrophobic_min_separation = 2;
  // 1/dist potential form
  //this->hydrophobic_min_cutoff = 2.0;
  //this->hydrophobic_max_cutoff = 8.0;
  //this->hydrophobic_max_Eshift = 0.125;
  // Spline potential form
  //this->hydrophobic_r = 6.0;
  //this->hydrophobic_half_delta = 2.0;
  /* electrostatics */
  this->recip_dielectric_param = 0.;
  this->debye_length_param = 0.;
  this->electrostatic_min_separation = 2;
  /* side chain hydrogen bond parameters */
  this->sidechain_hbond_strength_s2b = 0.;
  this->sidechain_hbond_strength_b2s = 0.;
  this->sidechain_hbond_strength_s2s = 0.;
  this->sidechain_hbond_cutoff = 0.;
  this->sidechain_hbond_decay_width = 0.;
  this->sidechain_hbond_min_separation = 2;
  this->sidechain_hbond_angle_cutoff = -1.;
  /* secondary radius of gyration */
  this->srgy_param = 0.0;
  this->srgy_offset = 0.0;
  this->hphobic_srgy_param = 0.0;
  this->hphobic_srgy_offset = 0.0;
  /*sbond */
  this->Sbond_strength = 0;
  this->Sbond_distance = 2.2;
  this->Sbond_cutoff = 0.2;
  this->Sbond_dihedral_cutoff = 0.35;

  /* fixed amino acids */
  this->fixed_aalist_file = NULL;

  /* external potential */
  this->external_potential_type = 0;
  for (int i=0; i<3; i++) {
    this->external_direction[i] = EXTERNAL_NONE;
    this->external_k[i] = 0.0;
    this->external_r0[i] = 0.0;
  }
  this->external_ztip = 0.0;
  this->external_constrained_aalist_file = NULL;
  this->external_potential_type2 = 0;
  for (int i=0; i<3; i++) {
    this->external_direction2[i] = EXTERNAL_NONE;
    this->external_k2[i] = 0.0;
    this->external_r02[i] = 0.0;
  }
  this->external_ztip2 = 0.0;
  this->external_constrained_aalist_file2 = NULL;

  //CAUTION!: aadict.c depends on params.c's model_params.  This means that
  //    initialize_sidechain_properties will have to be called after all updates
  //    of the vdW parameters; it can't be called from here, due to circular dependencies.
  this->sidechain_properties = calloc( 31, sizeof(sidechain_properties_) );
  /* vdw parameters might have changed */
  //initialize_sidechain_properties(this);

}

/* finalise model parameters */
void model_param_finalise(model_params *this) {

  /* gamma atoms */
  this->use_gamma_atoms = 0;
  this->use_original_gamma_atoms = 0;
  this->use_3_states = 0;
  this->fix_chi_angles = 0;
  this->fix_CA_atoms = 0;
  /* atomic radii */
  this->rca = 0.;
  this->rcb = 0.;
  this->rc = 0.;
  this->rn = 0.;
  this->ro = 0.;
  this->rs = 0.;
  this->rring = 0.;
  this->vdw_depth_ca = 0.;
  this->vdw_depth_cb = 0.;
  this->vdw_depth_c = 0.;
  this->vdw_depth_n = 0.;
  this->vdw_depth_o = 0.;
  this->vdw_depth_s = 0.;
  this->vdw_depth_ring = 0.;
  vdw_param_zero(this);
  if (this->vdw_gamma_gamma_cutoff) free(this->vdw_gamma_gamma_cutoff);
  if (this->vdw_gamma_nongamma_cutoff) free(this->vdw_gamma_nongamma_cutoff);
  this->vdw_backbone_cutoff = 0;
  this->vdw_use_extended_cutoff = 0;
  this->vdw_extended_cutoff = 0;
  this->vdw_potential = 0;
  this->vdw_function = NULL;
  this->vdw_clash_energy_at_hard_cutoff = 0;
  this->vdw_lj_neighbour_hard = 0;
  this->vdw_lj_hbonded_hard = 0;

  /* stress */
  this->stress_k = 0.;
  this->stress_angle = 0.;

  /* hydrogen bond */
  this->hboh2 = 0.;
  this->hbohn = 0.;
  this->hbcoh = 0.;
  //this->hboh_decay_width = 0.;
  //this->hbohn_decay_width = 0.;
  //this->hbcoh_decay_width = 0.;
  this->hbs = 0.;
  /* contact parameters */
  this->touch2 = 0.;
  this->part = 0.;
  this->split = 0.;
  this->sts = 0.;
  /* biasing force constants */
  if (this->contact_map_file) free(this->contact_map_file);
  this->bias_eta_beta = 0.;
  this->bias_eta_alpha = 0.;
  this->bias_kappa_alpha_3 = 0.;
  this->bias_kappa_alpha_4 = 0.;
  this->bias_kappa_beta = 0.;
  this->prt = 0.;
  this->bias_r_alpha = 0.;
  this->bias_r_beta = 0.;
  /* hydrophobicity */
  this->kauzmann_param = 0.;
  this->hydrophobic_cutoff_range = 0.;
  this->hydrophobic_min_separation = 0;
  //this->hydrophobic_min_cutoff = 0.;
  //this->hydrophobic_max_cutoff = 0.;
  //this->hydrophobic_max_Eshift = 0.;
  //this->hydrophobic_r = 0.;
  //this->hydrophobic_half_delta = 0.;
  /* electrostatics */
  this->recip_dielectric_param = 0.;
  this->debye_length_param = 0.;
  this->electrostatic_min_separation = 0;
  /* side chain hydrogen bond parameters */
  this->sidechain_hbond_strength_s2b = 0.;
  this->sidechain_hbond_strength_b2s = 0.;
  this->sidechain_hbond_strength_s2s = 0.;
  this->sidechain_hbond_cutoff = 0.;
  this->sidechain_hbond_decay_width = 0.;
  this->sidechain_hbond_min_separation = 0;
  this->sidechain_hbond_angle_cutoff = -0.;
  /* secondary radius of gyration */
  this->srgy_param = 0.;
  this->srgy_offset = 0.;
  this->hphobic_srgy_param = 0.;
  this->hphobic_srgy_offset = 0.;

  /* fixed amino acids */
  if (this->fixed_aalist_file) free(this->fixed_aalist_file);

  /* external potential */
  this->external_potential_type = 0;
  for (int i=0; i<3; i++) {
    this->external_direction[i] = EXTERNAL_NONE;
    this->external_k[i] = 0.0;
    this->external_r0[i] = 0.0;
  }
  this->external_ztip = 0.0;
  if (this->external_constrained_aalist_file) free(this->external_constrained_aalist_file);
  /* external potential #2 */
  this->external_potential_type2 = 0;
  for (int i=0; i<3; i++) {
    this->external_direction2[i] = EXTERNAL_NONE;
    this->external_k2[i] = 0.0;
    this->external_r02[i] = 0.0;
  }
  this->external_ztip2 = 0.0;
  if (this->external_constrained_aalist_file2) free(this->external_constrained_aalist_file2);

  if (this->sidechain_properties) free(this->sidechain_properties);
}


void flex_params_initialise(FLEX_params *this){
  this->number_of_processors = 0;
  this->flex_cmd = NULL;
  this->output_path = NULL;
  this->outputpdb_filename = NULL;
  this->output_path = realloc(this->output_path,DEFAULT_SHORT_STRING_LENGTH);
  strcpy(this->output_path,"/tmp/");
  this->outputpdb_filename = realloc(this->outputpdb_filename,DEFAULT_SHORT_STRING_LENGTH);
  strcpy(this->outputpdb_filename, "out.pdb");
  this->fromnode = 7;
  this->tonode = 7;
  this->step = 0.01;
  this->freq = 100;
  this->totalconf = 100;
  this->flex_dir = NULL;
  this->filenames_to_read_in = NULL;
  this->size_of_filename_to_read_in = 0;
  this->hstrength_cutoff = 1.0;
  this->only_bias_hbonds = 0;

  this->acceptance_rate_aim = 0.5;
  this->acceptance_rate_tolerance = 0.03;
  this->flex_changing_factor = 0.8;

  this->MCiter = 5000;

}

void flex_params_finalise(FLEX_params *this){
  this->number_of_processors = 0;
  this->fromnode = 7;
  this->tonode = 7;
  if(this->output_path){
    free(this->output_path);
    this->output_path = NULL;
  }
  if(this->outputpdb_filename){
    free (this->outputpdb_filename);
    this->outputpdb_filename = NULL;
  }
  if(this->flex_cmd){
    free (this->flex_cmd);
    this->flex_cmd = NULL;
  }
  if(this->flex_dir){
    free (this->flex_dir);
    this->flex_dir = NULL;
  }
  this->step = 0.01;
  this->freq = 100;
  this->totalconf = 100;

  if(this->size_of_filename_to_read_in != 0){
    int i;
    for(i = 0; i < this->size_of_filename_to_read_in;i++){
      if(this->filenames_to_read_in[i]) {
        free(this->filenames_to_read_in[i]);
        this->filenames_to_read_in[i] = NULL;
      }
    }
    free(this->filenames_to_read_in);
    this->filenames_to_read_in = NULL;
    this->size_of_filename_to_read_in = 0;

  }

  this->hstrength_cutoff = 1.0;
  this->only_bias_hbonds = 0;

  this->acceptance_rate_aim = 0.1;
  this->acceptance_rate_tolerance = 0.03;
  this->flex_changing_factor = 0.8;

}



/* vdW depths and shifts at the cutoffs for the atomic pairs */
void vdw_param_zero(model_params *this) {

   /* depths */

   this->vdw_depth_ca_ca = 0;
   this->vdw_depth_cb_cb = 0;
   this->vdw_depth_c_c   = 0;
   this->vdw_depth_n_n   = 0;
   this->vdw_depth_o_o   = 0;

   this->vdw_depth_ca_cb = 0;
   this->vdw_depth_ca_c  = 0;
   this->vdw_depth_ca_n  = 0;
   this->vdw_depth_ca_o  = 0;

   this->vdw_depth_cb_c  = 0;
   this->vdw_depth_cb_n  = 0;
   this->vdw_depth_cb_o  = 0;

   this->vdw_depth_c_n   = 0;
   this->vdw_depth_c_o   = 0;

   this->vdw_depth_n_o   = 0;

   /* shifts */
   this->vdw_shift = 0;

   this->vdw_Eshift_ca_ca = 0;
   this->vdw_Eshift_cb_cb = 0;
   this->vdw_Eshift_c_c   = 0;
   this->vdw_Eshift_n_n   = 0;
   this->vdw_Eshift_o_o   = 0;

   this->vdw_Eshift_ca_cb = 0;
   this->vdw_Eshift_ca_c  = 0;
   this->vdw_Eshift_ca_n  = 0;
   this->vdw_Eshift_ca_o  = 0;

   this->vdw_Eshift_cb_c  = 0;
   this->vdw_Eshift_cb_n  = 0;
   this->vdw_Eshift_cb_o  = 0;

   this->vdw_Eshift_c_n   = 0;
   this->vdw_Eshift_c_o   = 0;

   this->vdw_Eshift_n_o   = 0;

}


/***********************************************************/
/****                DEFAULT  PARAMETERS                ****/
/***********************************************************/


/* Set CD learnt values as default for all parameters learnt
   when using the LJ vdW interaction potential */
/* TODO: update parameters */
void set_lj_default_params(model_params *this) {

  /* The vdW cutoff distances have to be calculated later, because
     the cutoff calculating routine depends on vdw.c */
  this->vdw_potential = LJ_VDW_POTENTIAL;

  /* vdW  */
  /* atomic radii */
  /* radii from Ward et al, 1999 */
  // rca = 1.75, rcb = 1.75, rc = 1.65, rn = 1.55, ro = 1.40; 
  /* radii from Hopfinger, 1973 */
  // rca = 1.57, rcb = 1.57, rc = 1.42, rn = 1.29, ro = 1.29
  // rs = 1.8 from wikipedia I'm sorry ,NB
  /* LINUS */
  // rca = 1.85, rcb = 2.0, rc = 1.85, rn = 1.75, ro = 1.6, rs = 2.0; 
  this->rca = 2.43;
  this->rcb = 1.97;
  this->rc = 1.82;
  this->rn = 1.74;
  this->ro = 1.98;
  this->rs = 3.10;
  this->vdw_depth_ca = 0.018;
  this->vdw_depth_cb = 0.018;
  this->vdw_depth_c = 0.018;
  this->vdw_depth_n = 0.018;
  this->vdw_depth_o = 0.018;
  this->vdw_depth_s = 0.018;

  vdw_param_zero(this);
  vdw_param_calculate(this);

  this->vdw_clash_energy_at_hard_cutoff = 30; //default value for LJ

  /* stress */
  this->stress_k = 98.;

  /* hydrogen bond */
  this->hboh2 = 4.04;
  this->hbohn = 0.928;
  this->hbcoh = 0.772;
  this->hbs = 4.98;

  /* biasing force constants */
  this->bias_eta_beta = 3.7;
  this->bias_eta_alpha = 15.3;
  this->bias_kappa_beta = 0.85;
  this->bias_r_beta = 5.39;

  /* hydrophobicity */
  this->kauzmann_param = 0.122;
  //this->hydrophobic_cutoff_range = 2.8;
  //this->hydrophobic_min_separation = 2;

  /* electrostatics */
  /* side chain hydrogen bond parameters */
  /* secondary radius of gyration */

}


/* Set CD learnt values as default for all parameters learnt
   when using the hard cutoff vdW interaction potential */
/* TODO: update parameters */
void set_hard_cutoff_default_params(model_params *this) {

  /* The vdW cutoff distances have to be calculated later, because
     the cutoff calculating routine depends on vdw.c */
  this->vdw_potential = HARD_CUTOFF_VDW_POTENTIAL;

  /* vdW  */
  /* atomic radii */
  /* radii from Ward et al, 1999 */
  // rca = 1.75, rcb = 1.75, rc = 1.65, rn = 1.55, ro = 1.40; 
  /* radii from Hopfinger, 1973 */
  // rca = 1.57, rcb = 1.57, rc = 1.42, rn = 1.29, ro = 1.29
  // rs = 1.8 from wikipedia I'm sorry ,NB
  /* LINUS */
  // rca = 1.85, rcb = 2.0, rc = 1.85, rn = 1.75, ro = 1.6, rs = 2.0; 
  this->rca = 1.57;
  this->rcb = 1.57;
  this->rc = 1.42;
  this->rn = 1.29;
  this->ro = 1.29;
  this->rs = 2.00;
  this->vdw_depth_ca = 0;
  this->vdw_depth_cb = 0;
  this->vdw_depth_c = 0;
  this->vdw_depth_n = 0;
  this->vdw_depth_o = 0;
  this->vdw_depth_s = 0;

  vdw_param_zero(this);
  vdw_param_calculate(this);

  this->vdw_clash_energy_at_hard_cutoff = 10.4; //default value for LJ

  /* stress */
  this->stress_k = 90.;

  /* hydrogen bond */
  this->hboh2 = 4.05;
  this->hbohn = 0.930;
  this->hbcoh = 0.770;
  this->hbs = 4.95;

  /* biasing force constants */
  this->bias_eta_beta = 4.5;
  this->bias_eta_alpha = 18;
  this->bias_kappa_beta = 0.8;
  this->prt = 1.;
  this->bias_r_beta = 5.65;

  /* hydrophobicity */
  this->kauzmann_param = 0.13;
  //this->hydrophobic_cutoff_range = 2.8;
  //this->hydrophobic_min_separation = 2;

  /* electrostatics */
  /* side chain hydrogen bond parameters */
  /* secondary radius of gyration */

}


/***********************************************************/
/****                 COPYING  ROUTINES                 ****/
/***********************************************************/


/* Copy a string, making sure the null character is copied, too */
/* This will allocate new memory for the new string, so it has to have been freed before if it was used. */
void copy_string(char** const to, const char* const from){
  if (from==NULL) {
    *to=NULL;
  }
  else {
    //if (*to) free(*to); //this would crash things in model_params_copy
    //fprintf(stderr,"string from: %s, of length %d.\n",from,(int)strlen(from));
    *to = malloc((strlen(from)+1)*sizeof(char));
    if (*to == NULL) stop("Unable to allocate memory in copy_string.");
    strcpy(*to, from);
    //fprintf(stderr,"string from: %s, to: %s.\n",from,*to);
  }
}

/* Copy model parameters */
void model_params_copy(model_params *to, model_params *from) {
  *to = *from;
  /* strings */
  copy_string(&(to->contact_map_file), from->contact_map_file);
  copy_string(&(to->fixed_aalist_file), from->fixed_aalist_file);
  copy_string(&(to->external_constrained_aalist_file), from->external_constrained_aalist_file);
  copy_string(&(to->external_constrained_aalist_file2), from->external_constrained_aalist_file2);
  /* sidechain_properties */
  sidechain_properties_ *temp;
  temp = malloc(sizeof(sidechain_properties_) * 31);
  if (!temp) stop("Unable to allocate temp memory in model_params_copy.");
  to->sidechain_properties = temp;
  memcpy(to->sidechain_properties, from->sidechain_properties,
	31 * sizeof(sidechain_properties_));
  /* vdw cutoff matrices */
  if (from->vdw_gamma_gamma_cutoff) {
    double *temp1;
    temp1 = malloc(sizeof(double) * 702);
    if (!temp1) stop("Unable to allocate temp1 memory in model_params_copy.");
    to->vdw_gamma_gamma_cutoff = temp1;
    memcpy(to->vdw_gamma_gamma_cutoff, from->vdw_gamma_gamma_cutoff,
	702 * sizeof(double));
  }
  if (from->vdw_gamma_nongamma_cutoff) {
    double *temp2;
    temp2 = malloc(sizeof(double) * 702);
    if (!temp2) stop("Unable to allocate temp1 memory in model_params_copy.");
    to->vdw_gamma_nongamma_cutoff = temp2;
    memcpy(to->vdw_gamma_nongamma_cutoff, from->vdw_gamma_nongamma_cutoff,
	702 * sizeof(double));
  }
}


/* Copying simulation parameters (including model parameters) */
void sim_params_copy(simulation_params *to, simulation_params *from) {

  *to = *from;

  //TODO file pointers are copied, and files will be closed if the copy is finalised
  //to->infile = from->infile;
  //to->outfile = from->outfile;
  //to->checkpoint_file = from->checkpoint_file;

  //strings
  copy_string(&(to->prm), from->prm);
  copy_string(&(to->seq), from->seq);
  copy_string(&(to->sequence), from->sequence);
  copy_string(&(to->infile_name), from->infile_name);
  copy_string(&(to->outfile_name), from->outfile_name);
  copy_string(&(to->checkpoint_filename), from->checkpoint_filename);

  //double arrays
  double* temp = 0;
  temp = malloc(sizeof(double) * 36);
  if (!temp) stop("Unable to allocate temp memory in sim_params_copy.");
  to->energy_gradient = temp;
  memcpy(to->energy_gradient, from->energy_gradient, sizeof(double) * 36);

  temp = malloc(sizeof(double) * 36);
  if (!temp) stop("Unable to allocate temp memory in sim_params_copy.");
  to->energy_probe_1_this = temp;
  memcpy(to->energy_probe_1_this, from->energy_probe_1_this, sizeof(double) * 36);

  temp = malloc(sizeof(double) * 36);
  if (!temp) stop("Unable to allocate temp memory in sim_params_copy.");
  to->energy_probe_1_last = temp;
  memcpy(to->energy_probe_1_last, from->energy_probe_1_last, sizeof(double) * 36);

  //integer arrays
  int* temp2 = malloc(sizeof(int) * 36);
  if (!temp2) stop("Unable to allocate temp memory2 in sim_params_copy.");
  to->energy_probe_1_calc = temp2;
  memcpy(to->energy_probe_1_calc, from->energy_probe_1_calc, sizeof(int) * 36);

  if (from->seq != NULL && from->sequence != NULL && from->MC_lookup_table != NULL && from->MC_lookup_table_n != NULL) {
    int N = 4 * (from->NAA + strlen(from->sequence) - strlen(from->seq));
    temp2 = malloc(sizeof(int) * N);
    if (!temp2) stop("Unable to allocate temp memory2 in sim_params_copy.");
    to->MC_lookup_table = temp2;
    memcpy(to->MC_lookup_table, from->MC_lookup_table, sizeof(int) * N);
    //MC lookup table length
    temp2 = malloc(sizeof(int) * 4);
    if (!temp2) stop("Unable to allocate temp memory2 in sim_params_copy.");
    to->MC_lookup_table_n = temp2;
    memcpy(to->MC_lookup_table_n, from->MC_lookup_table_n, sizeof(int) * 4);
  } else {
    to->MC_lookup_table = NULL;
    to->MC_lookup_table_n = NULL;
  }

  //to->protein_model = malloc(sizeof(model_params));
  //if (!to->protein_model) stop("Unable to allocate protein_model memory in sim_params_copy."); 




  model_params_copy(&(to->protein_model), &(from->protein_model));
  flex_params_copy(&(to->flex_params), &(from->flex_params));
}


void flex_params_copy(FLEX_params *to, FLEX_params *from){
  *to = *from;
  copy_string(&(to->output_path), from->output_path);
  copy_string(&(to->outputpdb_filename), from->outputpdb_filename);
  copy_string(&(to->flex_cmd), from->flex_cmd);
  copy_string(&(to->flex_dir), from->flex_dir);

  if(to->size_of_filename_to_read_in != 0){
    to->filenames_to_read_in = NULL;
    to->filenames_to_read_in = (char**)malloc(sizeof(char*)*to->size_of_filename_to_read_in);
    int i;
    for(i = 0; i < to->size_of_filename_to_read_in; i++){
      copy_string(&(to->filenames_to_read_in[i]),from->filenames_to_read_in[i]);
    }

  }


}


/***********************************************************/
/****                 READING  ROUTINES                 ****/
/***********************************************************/

void flex_setup_command(FLEX_params *this){
  if(this->size_of_filename_to_read_in > 0){
    int f;
    for(f = 0; f <this->size_of_filename_to_read_in; f++){
      free(this->filenames_to_read_in[f]);
    }
    free(this->filenames_to_read_in);
    this->filenames_to_read_in = NULL;
  }


  this->size_of_filename_to_read_in = 0;
  this->mode = this->fromnode + (int)(rand()/(double)RAND_MAX * (this->tonode-this->fromnode+1)) % (this->tonode-this->fromnode+1);
  //for(mode = this->fromnode; mode <= this->tonode; mode++){
    this->freq = this->totalconf/2 + (int)(rand()/(double)RAND_MAX * (this->totalconf/2)) % (this->totalconf/2);
    //for(freq = this->freq; freq <= this->totalconf; freq+=this->freq){
      //this->freq

      this->size_of_filename_to_read_in+=2;
      this->filenames_to_read_in = (char**)realloc(this->filenames_to_read_in,sizeof(char*)*this->size_of_filename_to_read_in);
      this->filenames_to_read_in[this->size_of_filename_to_read_in-1] = NULL;
      this->filenames_to_read_in[this->size_of_filename_to_read_in-2] = NULL;
      this->filenames_to_read_in[this->size_of_filename_to_read_in-1] = (char*)malloc(sizeof(char)*DEFAULT_LONG_STRING_LENGTH);
      this->filenames_to_read_in[this->size_of_filename_to_read_in-2] = (char*)malloc(sizeof(char)*DEFAULT_LONG_STRING_LENGTH);
      sprintf(this->filenames_to_read_in[this->size_of_filename_to_read_in-2],"%sRuns/Mode%02d-pos/out_froda_%08d.pdb",this->output_path,this->mode,this->freq);
      sprintf(this->filenames_to_read_in[this->size_of_filename_to_read_in-1],"%sRuns/Mode%02d-neg/out_froda_%08d.pdb",this->output_path,this->mode,this->freq);

    //}

  //}

  //fprintf(stdout,"Mode: %d Freq %d\n",mode,freq); fflush(stdout);

  if(this->flex_cmd){free(this->flex_cmd); this->flex_cmd = NULL;}
  this->flex_cmd = (char*)malloc(sizeof(char)*1000);
  sprintf(this->flex_cmd,"%s%s %s %s %d %d %lf %d %d %s",this->flex_dir,"flex_script.sh ",this->output_path,this->outputpdb_filename,this->mode,this->mode,this->step, this->freq,this->freq,this->flex_dir);

}




/* read in model parameters */
int flex_param_read(char *prm, FLEX_params *this){
  char outdir_string[256]="";
  char flex_dir[256]="";
  int k = sscanf(prm,"FLEX=%d,%d,%d,%d,%256[^,],%256[^,],%lf,%lf,%lf",&(this->number_of_processors),&(this->totalconf),&(this->fromnode),&(this->tonode),outdir_string,flex_dir,&(this->acceptance_rate_aim),&(this->hstrength_cutoff),&(this->step));

  int returnval = 5;

  //if(this->number_of_processors == 0) return 4;
  if(this->hstrength_cutoff > 1.0 || this->hstrength_cutoff <= 0.0){
    stop("FLEX hbond strength cutoff must be in (0,1]\n");
  }
  if(this->fromnode < 7 || this->tonode < 7){
     stop("Error: node numbers for FLEX must be at least 7\n");
  }
  this->flex_dir = (char*)malloc(sizeof(char)*DEFAULT_SHORT_STRING_LENGTH);

  if(this->totalconf == 0) stop("Error totalconf = 0\n");

  if(k>=5){
    sprintf(this->output_path,"%s",outdir_string);
    returnval += strlen(outdir_string);

    for (int digit=1; digit<=this->totalconf; digit*=10) {
      returnval ++;
    }
    for (int digit=1; digit<=this->fromnode; digit*=10) {
      returnval ++;
    }
    for (int digit=1; digit<=this->tonode; digit*=10) {
      returnval ++;
    }
    for (int digit=1; digit<=this->number_of_processors; digit*=10) {
          returnval ++;
    }

    returnval += 4;
  }
  if(k >= 6){
    strcpy(this->flex_dir,flex_dir);
	returnval += strlen(flex_dir) + 1;
  }

  flex_setup_command(this);

  return returnval;
}

void model_param_read(char *prm, /* input line */
		model_params *this,FLEX_params *nma_params) {

    int k;
    int start;
    int explicit_contact_map_file = 0;
    int found_param;
    char error_string[DEFAULT_LONG_STRING_LENGTH]="";

    if (prm == NULL) {
	fprintf(stderr,"WARNING! No command line parameter line was given.\n");
	return;
    }
    //fprintf(stderr,"%s\n",prm);

    /* Read parameters one by one */
    while (strlen(prm)>0) {

	//fprintf(stderr,"%s\n",prm);

	start = 1;
	found_param = 0;

	/* backbone H-bond */
	k=sscanf(prm, "Hbond=%lf,%lf,%lf,%lf",
	//k=sscanf(prm, "Hbond=%lf,%lf,%lf,%lf,%lf,%lf,%lf",
		&(this->hboh2),
		&(this->hbohn),
		&(this->hbcoh),
		&(this->hbs) /*,
		&(this->hboh_decay_width),
		&(this->hbohn_decay_width),
		&(this->hbcoh_decay_width) */ );
	//fprintf(stderr,"Read in %d h-bond params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 6;
	}

	/* contact parameters */
	k=sscanf(prm, "Contact=%lf,%lf,%lf,%lf",
		&(this->touch2),
		&(this->part),
		&(this->split),
		&(this->sts));
	//fprintf(stderr,"Read in %d contact params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 8;
	}

	char gamma_string[256];
	/* gamma atoms of the model */
	k=sscanf(prm, "Gamma=%256[^,],%d,%d,%d,%d",
		gamma_string,
		&(this->use_original_gamma_atoms),
		&(this->use_3_states),
		&(this->fix_chi_angles),
		&(this->fix_CA_atoms));
	//fprintf(stderr,"Read in %d gamma params.\n",k);
	if (k>0) {
//		fprintf(stderr,"The gamma string is %s.\n",gamma_string);
		if (strcmp(gamma_string,"NONE")==0) {
			this->use_gamma_atoms = NO_GAMMA;
		} else if (strcmp(gamma_string,"LINUS_GAMMA")==0) {
			this->use_gamma_atoms = LINUS_GAMMA;
		} else if (strcmp(gamma_string,"CORRECT_GAMMA")==0) {
			this->use_gamma_atoms = CORRECT_GAMMA;
		} else if (strcmp(gamma_string,"CORRECT_KMQR_GAMMA")==0) {
			this->use_gamma_atoms = CORRECT_KMQR_GAMMA;
		} else {
			sprintf(error_string,"Unknown value for gamma model (%s).  It must be one of NONE, LINUS_GAMMA, CORRECT_GAMMA and CORRECT_KMQR_GAMMA.",gamma_string);
			stop(error_string);
		}
		if (this->use_gamma_atoms==NO_GAMMA) {
			fprintf(stderr,"WARNING! Not using gamma atoms.\n");
			fprintf(stderr,"WARNING! If you want to set an extended cutoff value for this model to make it faster, do not forget to specify if in the VDW=... parameters, otherwise the default value for the gamma atom including model will be used.\n");
		}

		found_param += 1;
		start = strlen(gamma_string)+6; //B=gamma_string...
	}

	/* vdW radii */
	char potential_string[256];
	k=sscanf(prm, "VDW=%256[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%lf,%lf",
		   potential_string,
		   &(this->rca),
		   &(this->rcb),
		   &(this->rc),
		   &(this->rn),
		   &(this->ro),
		   &(this->rs),
		   &(this->rel_vdw_cutoff),
		   &(this->vdw_depth_ca),
		   &(this->vdw_depth_cb),
		   &(this->vdw_depth_c),
		   &(this->vdw_depth_n),
		   &(this->vdw_depth_o),
		   &(this->vdw_depth_s),
		   &(this->vdw_uniform_depth),
		   &(this->vdw_use_extended_cutoff),
		   &(this->vdw_extended_cutoff),
		   &(this->vdw_clash_energy_at_hard_cutoff));
	if (k>0) {
		//fprintf(stderr,"%s\n",potential_string);
		if (strcmp(potential_string,"hard_cutoff")==0) {
			this->vdw_potential = HARD_CUTOFF_VDW_POTENTIAL;
			//set defaults for the hard cutoff model
			if (k<18) this->vdw_clash_energy_at_hard_cutoff = 10.4;
			if (k<14) {
				this->vdw_depth_ca = 0;
				this->vdw_depth_cb = 0;
				this->vdw_depth_c = 0;
				this->vdw_depth_n = 0;
				this->vdw_depth_o = 0;
				this->vdw_depth_s = 0;
			}
		} else if (strcmp(potential_string,"lj")==0) {
			this->vdw_potential = LJ_VDW_POTENTIAL;
		} else {
			sprintf(error_string,"Unknown value for vdW potential model (%s).  It must be one of hard_cutoff and lj.",potential_string);
			stop(error_string);
		}
		found_param += 1;
		start = strlen(potential_string) + 4;
	}

	vdw_param_calculate(this);
	//CAUTION!: aadict.c depends on params.c's model_params.  This means that
	//    initialize_sidechain_properties will have to be called after all updates
	//    of the vdW parameters; it can't be called from here, due to circular dependencies.
	//initialize_sidechain_properties(this);

	/* LJ parameters */
	//k=sscanf(prm, "LJ=%d,%d",
	//	   potenti,
	//	   &(this->vdw_lj_neighbour_hard),
	//	   &(this->vdw_lj_hbonded_hard));
	//if (k>0) {
	//	found_param += 1;
	//	start = 3;
	//}

	/* stress */
	k=sscanf(prm, "Stress=%lf,%lf",
		   &(this->stress_k),
		   &(this->stress_angle));
	if (k>0) {
		found_param += 1;
		start = 7;
	}

	/* bias potential */
	char contact_map_file[256];
	if ((k=sscanf(prm, "Bias=%256[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
		   contact_map_file,
		   &(this->bias_eta_beta),
		   &(this->bias_eta_alpha),
		   &(this->bias_kappa_alpha_3),
		   &(this->bias_kappa_alpha_4),
		   &(this->bias_kappa_beta),
		   &(this->prt),
		   &(this->bias_r_alpha),
		   &(this->bias_r_beta))) > 0) {
		this->contact_map_file = realloc(this->contact_map_file,DEFAULT_SHORT_STRING_LENGTH);
		strcpy(this->contact_map_file,contact_map_file);
//		this->contact_map_file = contact_map_file;
		explicit_contact_map_file = 1;

		found_param += 1;
		start = strlen(contact_map_file)+5; //B=contact_map_file...
	}
	//fprintf(stderr,"Read in %d bias params.\n",k);

	/* hydrophobic parameters */
	// 1/dist potential form
	//k=sscanf(prm, "Hydrophobic=%lf,%d,%lf",
	//	&(this->kauzmann_param),
	//	&(this->hydrophobic_min_separation),
	//	&(this->hydrophobic_max_cutoff));
	//if (k==3) this->hydrophobic_max_Eshift = 1.0 / this->hydrophobic_max_cutoff;
	// Spline potential form
	//k=sscanf(prm, "Hydrophobic=%lf,%d,%lf,%lf",
	//	&(this->kauzmann_param),
	//	&(this->hydrophobic_min_separation),
	//	&(this->hydrophobic_r),
	//	&(this->hydrophobic_half_delta));
	// Linear decay potential form (default)
	k=sscanf(prm, "Hydrophobic=%lf,%d,%lf",
		&(this->kauzmann_param),
		&(this->hydrophobic_min_separation),
		&(this->hydrophobic_cutoff_range));
	//fprintf(stderr,"Read in %d hydrophobic params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 12;
	}
	
	k=sscanf(prm, "Fixit=%d",
		&(this->fixit));
	if (k>0) {
		found_param += 1;
		start = 6;
	}

	/* electrostatic parameters */
	k=sscanf(prm, "Elec=%lf,%lf,%d",
		&(this->recip_dielectric_param),
		&(this->debye_length_param),
		&(this->electrostatic_min_separation));
	//fprintf(stderr,"Read in %d electrostatic params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 5;
	}

	/* side chain hydrogen bond parameters */
	k=sscanf(prm, "SchHbond=%lf,%lf,%lf,%lf,%lf,%d",
		&(this->sidechain_hbond_strength_s2b),
		&(this->sidechain_hbond_strength_b2s),
		&(this->sidechain_hbond_strength_s2s),
		&(this->sidechain_hbond_angle_cutoff),
		&(this->sidechain_hbond_decay_width),
		&(this->sidechain_hbond_min_separation));
	//fprintf(stderr,"Read in %d sidechain h-bond params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 9;
	}

	/* secondary radius of gyration parameters */
	k=sscanf(prm, "Rgyr=%lf,%lf,%lf,%lf",
		&(this->srgy_param),
		&(this->srgy_offset),
		&(this->hphobic_srgy_param),
		&(this->hphobic_srgy_offset));
	//fprintf(stderr,"Read in %d secondary radius of gyration params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 5;

	}

	/* S-S bond */
	k=sscanf(prm,"SSbond=%lf,%lf,%lf,%lf",
		&(this->Sbond_strength),
		&(this->Sbond_distance),
		&(this->Sbond_cutoff),
		&(this->Sbond_dihedral_cutoff));
	if (k>0) {
		found_param += 1;
		start = 7;
	}

	/* fixed amino acids */
	char fixed_file[256];
	if ((k=sscanf(prm, "fixed=%256[^,]",fixed_file)) == 1) {
		fprintf(stderr,"setting fixed amino acid list file %s\n",fixed_file);
		this->fixed_aalist_file = realloc(this->fixed_aalist_file,DEFAULT_SHORT_STRING_LENGTH);
		strcpy(this->fixed_aalist_file,fixed_file);
		found_param += 1;
		start = strlen(fixed_file) + strlen("fixed=");
	}


	/* external potential */
	int xdir = 0;
	int ydir = 0;
	int zdir = 0;
	char constraint_file[256];
	k=sscanf(prm,"external=%d",&(this->external_potential_type));
	//fprintf(stderr, "haha found external potential parameters\n");
	if (k>0) {
		fprintf(stderr,"found external potential parameters %d \n", k);
		int l;
		if (this->external_potential_type == 1) {
			l=sscanf(prm+9,"%d,%256[^,],%lf,%lf",
				&(this->external_potential_type), constraint_file, &(this->external_k[0]), &(this->external_r0[0])); //save r0, k into 1st vector element
			fprintf(stderr, "found external potential parameters %d \n", l);
			//stop("unimplemented");
			if (l>1) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file);
				this->external_constrained_aalist_file = realloc(this->external_constrained_aalist_file,DEFAULT_SHORT_STRING_LENGTH);
				strcpy(this->external_constrained_aalist_file,constraint_file);
			}

		} else if (this->external_potential_type == 2) { /* unimplemented */
			stop("unimplemented");
			l=sscanf(prm+9,"%d,%d,%lf,%lf,%d,%lf,%lf,%d,%lf,%lf",
				&(this->external_potential_type), &xdir, &(this->external_k[0]), &(this->external_r0[0]),
				&ydir, &(this->external_k[1]), &(this->external_r0[1]),
				&zdir, &(this->external_k[2]), &(this->external_r0[2]));
			if (l>0) {
				/* check x direction */
				if (xdir == -1) this->external_direction[0] = EXTERNAL_NEGATIVE;
				else if (xdir == 1) this->external_direction[0] = EXTERNAL_POSITIVE;
				else if (xdir == 2) this->external_direction[0] = EXTERNAL_POSNEG;
				else if (xdir == 0) this->external_direction[0] = EXTERNAL_NONE;
				else stop("external direction has to be one of -1(negative), 0(none), 1(positive) or 2 (positive and negative)");
				/* check y direction */
				if (ydir == -1) this->external_direction[1] = EXTERNAL_NEGATIVE;
				else if (ydir == 1) this->external_direction[1] = EXTERNAL_POSITIVE;
				else if (ydir == 2) this->external_direction[1] = EXTERNAL_POSNEG;
				else if (ydir == 0) this->external_direction[1] = EXTERNAL_NONE;
				else stop("external direction has to be one of -1(negative), 0(none), 1(positive) or 2 (positive and negative)");
				/* check z direction */
				if (zdir == -1) this->external_direction[2] = EXTERNAL_NEGATIVE;
				else if (zdir == 1) this->external_direction[2] = EXTERNAL_POSITIVE;
				else if (zdir == 2) this->external_direction[2] = EXTERNAL_POSNEG;
				else if (zdir == 0) this->external_direction[2] = EXTERNAL_NONE;
				else stop("external direction has to be one of -1(negative), 0(none), 1(positive) or 2 (positive and negative)");

			}
		} else if (this->external_potential_type == 3) { /* conical potential */
			l=sscanf(prm+9,"%d,%256[^,],%lf,%lf,%lf",
				&(this->external_potential_type), constraint_file, &(this->external_k[0]), &(this->external_r0[0]),&(this->external_ztip)); //save r0, k into 1st vector element
			if (l==5) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file);
				this->external_constrained_aalist_file = realloc(this->external_constrained_aalist_file,DEFAULT_SHORT_STRING_LENGTH);
				strcpy(this->external_constrained_aalist_file,constraint_file);
			} else {
				stop("conical potential needs 4 parameters (constraint list file,k,r0,ztip)");
			}
		}
		else if (this->external_potential_type == 5) {
			l = sscanf(prm + 9, "%d,%256[^,],%lf,%lf",
				&(this->external_potential_type), constraint_file, &(this->external_k[0]), &(this->external_r0[0])); //save r0, k into 1st vector element
			fprintf(stderr, "found external potential parameters %d \n", l);
			//stop("unimplemented");
			if (l > 1) {
				fprintf(stderr, "setting constrained list file %s\n", constraint_file);
				this->external_constrained_aalist_file = realloc(this->external_constrained_aalist_file, DEFAULT_SHORT_STRING_LENGTH);
				strcpy(this->external_constrained_aalist_file, constraint_file);
			}
		} else {
			stop("unknown external potential type, has to be 1 (cylindrical around z) or 3 (conical around z).");
		}
		found_param += 1;
		start = strlen(constraint_file) + 9 + 2; //2 comes from the type and the comma after
		//fprintf(stderr,"setting found_param to %d and start to %d\n",found_param,start);
	}
	/* external potential 2 */
	char constraint_file2[256];
	k=sscanf(prm,"external2=%d",&(this->external_potential_type2));
	if (k>0) {
		fprintf(stderr,"found external2 potential parameters %d\n", this->external_potential_type2);
		int l;
		if (this->external_potential_type2 == 1) {
			l=sscanf(prm+10,"%d,%256[^,],%lf,%lf",
				&(this->external_potential_type2), constraint_file2, &(this->external_k2[0]), &(this->external_r02[0])); //save r0, k into 1st vector element
			if (l>1) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file2);
				this->external_constrained_aalist_file2 = realloc(this->external_constrained_aalist_file2,DEFAULT_SHORT_STRING_LENGTH);
				strcpy(this->external_constrained_aalist_file2,constraint_file2);
			}
		} else if (this->external_potential_type2 == 3) { /* conical potential */
			l=sscanf(prm+9,"%d,%256[^,],%lf,%lf,%lf",
				&(this->external_potential_type2), constraint_file2, &(this->external_k2[0]), &(this->external_r02[0]),&(this->external_ztip2)); //save r0, k into 1st vector element
			if (l==5) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file2);
				this->external_constrained_aalist_file2 = realloc(this->external_constrained_aalist_file2,DEFAULT_SHORT_STRING_LENGTH);
				strcpy(this->external_constrained_aalist_file2,constraint_file2);
			} else {
				stop("conical potential needs 4 parameters (constraint list file,k,r0,ztip)");
			}
		} else if (this->external_potential_type2 == 4) {
			l = sscanf(prm + 10, "%d,%256[^,],%lf,%lf", &(this->external_potential_type2), constraint_file2, &(this->external_k2[0]), &(this->external_r02[0]));
			fprintf(stderr, "setting cyclic peptide %s \n", constraint_file2);
		} else {
			stop("unknown external2 potential type, has to be 1.");
		}
		found_param += 1;
		start = strlen(constraint_file2) + 10 + 2; //2 comes from the type and the comma after
		//fprintf(stderr,"setting found_param to %d and start to %d\n",found_param,start);
	}

    /*dummy position for FLEX */
	char temp[DEFAULT_LONG_STRING_LENGTH];
	k = sscanf(prm,"FLEX=%s",temp);
	if(k>0){
	  start = flex_param_read(prm, nma_params);
	  found_param += 1;

	}



	/* Check for more than one parameters found */
	if (found_param > 1) {
		sprintf(error_string,"The same section of the parameter line was used to read in more than one parameters.  Prm: %s",prm);
		stop(error_string);
	}

	/* Check for unknown parameters */
	if (found_param == 0 ) {
		sprintf(error_string,"Unknown parameter in prm string: %s",prm);
		stop(error_string);
	}

	/* Shift the prm to the next parameter entry */
	int next = 0;
	for (int i=start; (next==0 && i<=strlen(prm)); i++) {
		if (isalpha(prm[i])) next = i;
	}
	if (next > 0 && next < strlen(prm)) {
		/* the new prm is >=1 long */
		prm = prm+next;
	} else {
		//model_param_print(*this);
		if (!explicit_contact_map_file) fprintf(stderr,"WARNING: No contact map specified.");
		return;
	}

    }

}




/* calculating vdW depths and shifts at the cutoffs for the atomic pairs,
   not to have to recalculate on the fly */
void vdw_param_calculate(model_params *this) {

   /* depths */

   // quick check that nothing is negative
   if (this->vdw_depth_ca < 0) this->vdw_depth_ca = 0;
   if (this->vdw_depth_cb < 0) this->vdw_depth_cb = 0;
   if (this->vdw_depth_c  < 0) this->vdw_depth_c  = 0;
   if (this->vdw_depth_n  < 0) this->vdw_depth_n  = 0;
   if (this->vdw_depth_o  < 0) this->vdw_depth_o  = 0;
   if (this->vdw_depth_s  < 0) this->vdw_depth_s  = 0;

   if (this->vdw_uniform_depth) {
   this->vdw_depth_ca_sqrt = this->vdw_depth_cb_sqrt =
			     this->vdw_depth_c_sqrt  =
			     this->vdw_depth_n_sqrt  =
			     this->vdw_depth_o_sqrt  =
			     this->vdw_depth_s_sqrt  = sqrt(this->vdw_depth_ca);
   this->vdw_depth_ca_ca = this->vdw_depth_cb_cb =
			   this->vdw_depth_c_c   =
			   this->vdw_depth_n_n   =
			   this->vdw_depth_o_o   =
			   this->vdw_depth_ca_cb =
			   this->vdw_depth_ca_c  =
			   this->vdw_depth_ca_n  =
			   this->vdw_depth_ca_o  =

			   this->vdw_depth_cb_c  =
			   this->vdw_depth_cb_n  =
			   this->vdw_depth_cb_o  =

			   this->vdw_depth_c_n   =
			   this->vdw_depth_c_o   =

			   this->vdw_depth_n_o   = this->vdw_depth_ca;
   /* shifts */

   double shift = pow(1/this->rel_vdw_cutoff,12) - 2 * pow(1/this->rel_vdw_cutoff,6);
   this->vdw_shift = shift;

   this->vdw_Eshift_ca_ca = this->vdw_Eshift_cb_cb =
			    this->vdw_Eshift_c_c   =
			    this->vdw_Eshift_n_n   =
			    this->vdw_Eshift_o_o   =
			    this->vdw_Eshift_ca_cb =
			    this->vdw_Eshift_ca_c  =
			    this->vdw_Eshift_ca_n  =
			    this->vdw_Eshift_ca_o  =
			    this->vdw_Eshift_cb_c  =
			    this->vdw_Eshift_cb_n  =
			    this->vdw_Eshift_cb_o  =
			    this->vdw_Eshift_c_n   =
			    this->vdw_Eshift_c_o   =
			    this->vdw_Eshift_n_o   = this->vdw_depth_ca_ca * shift;
   } else {
   this->vdw_depth_ca_sqrt = sqrt(this->vdw_depth_ca);
   this->vdw_depth_cb_sqrt = sqrt(this->vdw_depth_cb);
   this->vdw_depth_c_sqrt = sqrt(this->vdw_depth_c);
   this->vdw_depth_n_sqrt = sqrt(this->vdw_depth_n);
   this->vdw_depth_o_sqrt = sqrt(this->vdw_depth_o);
   this->vdw_depth_s_sqrt = sqrt(this->vdw_depth_s);
  
   this->vdw_depth_ca_ca = this->vdw_depth_ca;
   this->vdw_depth_cb_cb = this->vdw_depth_cb;
   this->vdw_depth_c_c   = this->vdw_depth_c;
   this->vdw_depth_n_n   = this->vdw_depth_n;
   this->vdw_depth_o_o   = this->vdw_depth_o;

   this->vdw_depth_ca_cb = sqrt(this->vdw_depth_ca * this->vdw_depth_cb);
   this->vdw_depth_ca_c  = sqrt(this->vdw_depth_ca * this->vdw_depth_c );
   this->vdw_depth_ca_n  = sqrt(this->vdw_depth_ca * this->vdw_depth_n );
   this->vdw_depth_ca_o  = sqrt(this->vdw_depth_ca * this->vdw_depth_o );

   this->vdw_depth_cb_c  = sqrt(this->vdw_depth_cb * this->vdw_depth_c );
   this->vdw_depth_cb_n  = sqrt(this->vdw_depth_cb * this->vdw_depth_n );
   this->vdw_depth_cb_o  = sqrt(this->vdw_depth_cb * this->vdw_depth_o );

   this->vdw_depth_c_n   = sqrt(this->vdw_depth_c  * this->vdw_depth_n );
   this->vdw_depth_c_o   = sqrt(this->vdw_depth_c  * this->vdw_depth_o );

   this->vdw_depth_n_o   = sqrt(this->vdw_depth_n  * this->vdw_depth_o );

   /* shifts */

   double shift = pow(1/this->rel_vdw_cutoff,12) - 2 * pow(1/this->rel_vdw_cutoff,6);
   this->vdw_shift = shift;

   this->vdw_Eshift_ca_ca = this->vdw_depth_ca_ca * shift;
   this->vdw_Eshift_cb_cb = this->vdw_depth_cb_cb * shift;
   this->vdw_Eshift_c_c   = this->vdw_depth_c_c   * shift;
   this->vdw_Eshift_n_n   = this->vdw_depth_n_n   * shift;
   this->vdw_Eshift_o_o   = this->vdw_depth_o_o   * shift;

   this->vdw_Eshift_ca_cb = this->vdw_depth_ca_cb * shift;
   this->vdw_Eshift_ca_c  = this->vdw_depth_ca_c  * shift;
   this->vdw_Eshift_ca_n  = this->vdw_depth_ca_n  * shift;
   this->vdw_Eshift_ca_o  = this->vdw_depth_ca_o  * shift;

   this->vdw_Eshift_cb_c  = this->vdw_depth_cb_c  * shift;
   this->vdw_Eshift_cb_n  = this->vdw_depth_cb_n  * shift;
   this->vdw_Eshift_cb_o  = this->vdw_depth_cb_o  * shift;

   this->vdw_Eshift_c_n   = this->vdw_depth_c_n   * shift;
   this->vdw_Eshift_c_o   = this->vdw_depth_c_o   * shift;

   this->vdw_Eshift_n_o   = this->vdw_depth_n_o   * shift;

  }

}


/***********************************************************/
/****                 PRINTING ROUTINES                 ****/
/***********************************************************/


/* Printing the vdW backbone cutoff and the gamma - gamma and gamma - nongamma cutoff matrices */
void print_vdw_cutoff_distances(model_params *mod_params, FILE *outfile) {

	int l, m;

	fprintf(outfile,"vdW backbone cutoff: %g\n",mod_params->vdw_backbone_cutoff);

	if (mod_params->use_gamma_atoms == NO_GAMMA) {
	    fprintf(outfile,"vdW cutoff matrices are not used (no gamma atoms in the model).\n");
	    return;
	}

	if (mod_params->vdw_gamma_gamma_cutoff == NULL) {
	    fprintf(outfile,"vdW gamma - gamma cutoff was not allocated.\n");
	    return;
	}
	if (mod_params->vdw_gamma_nongamma_cutoff == NULL) {
	    fprintf(outfile,"vdW gamma - gamma cutoff was not allocated.\n");
	    return;
	}

	/* Print the cutoff matrices */
	/* The 26x26 array of gamma - gamma contact cutoffs */
	fprintf(outfile,"const double maxvdw_gamma_gamma[702] ={ ");
	for(l = 0; l < 27; l++){
	    for(m = 0; m < 26; m++){
		fprintf(outfile,"%f, ",mod_params->vdw_gamma_gamma_cutoff[l*26+m]);
	    }
	    fprintf(outfile,"\n");
	}
	fprintf(outfile," };\n");

	/* The 26x26 array of gamma - nongamma contact cutoffs */
	fprintf(outfile,"const double maxvdw_gamma_nongamma[702] ={ ");
	for(l = 0; l < 27; l++){
	    for(m = 0; m < 26; m++){
		fprintf(outfile,"%f, ",mod_params->vdw_gamma_nongamma_cutoff[l*26+m]);
	    }
	    fprintf(outfile,"\n");
	}
	fprintf(outfile," };\n");

}

/* print model parameters */
void model_param_print(model_params this, FILE *outfile) {

  fprintf(outfile,"================MODEL=PARAMETERS===================\n");
  /* gamma atoms */
  fprintf(outfile,"use gamma atoms?(no:%d,LINUS:%d,correct:%d,LINUS+correct_KMQR:%d) %d\n",NO_GAMMA,LINUS_GAMMA,CORRECT_GAMMA,CORRECT_KMQR_GAMMA,this.use_gamma_atoms);
  fprintf(outfile,"use original sidechain dihedrals when reading (1), or pick a random one(0)? %d\n",this.use_original_gamma_atoms);
  fprintf(outfile,"adjust to 3 states when reading in? %d\n",this.use_3_states);
  fprintf(outfile,"fix the side chain dihedral angles in the mc? %d\n",this.fix_chi_angles);
  fprintf(outfile,"fix CA atoms in the mc (only pivot peptide bonds)? %d\n",this.fix_CA_atoms);
  /* atomic radii */
  fprintf(outfile,"VDW RADII\n");
  fprintf(outfile,"r(CA_) %g\n",this.rca);
  fprintf(outfile,"r(CB_) %g\n",this.rcb);
  fprintf(outfile,"r(C__) %g\n",this.rc);
  fprintf(outfile,"r(O__) %g\n",this.ro);
  fprintf(outfile,"r(N__) %g\n",this.rn);
  fprintf(outfile,"r(S__) %g\n",this.rs);
  fprintf(outfile,"rel. cutoff %g\n",this.rel_vdw_cutoff);
  fprintf(outfile,"uniform depth? %d\n",this.vdw_uniform_depth);
  fprintf(outfile,"eps(CA_) %g\n",this.vdw_depth_ca);
  fprintf(outfile,"eps(CB_) %g\n",this.vdw_depth_cb);
  fprintf(outfile,"eps(C__) %g\n",this.vdw_depth_c);
  fprintf(outfile,"eps(N__) %g\n",this.vdw_depth_n);
  fprintf(outfile,"eps(O__) %g\n",this.vdw_depth_o);
  fprintf(outfile,"eps(S__) %g\n",this.vdw_depth_s);
  fprintf(outfile,"vdW potential model (%d: hard cutoff, %d: LJ) %d\n",HARD_CUTOFF_VDW_POTENTIAL,LJ_VDW_POTENTIAL,this.vdw_potential);
  fprintf(outfile,"Use hard cutoff part of LJ for neighbouring amino acid pairs? (0:no, 1:yes) %d\n",this.vdw_lj_neighbour_hard);
  fprintf(outfile,"Use hard cutoff part of LJ for H-bonded amino acid pairs (and their neighbours)? (0:no, 1:yes) %d\n",this.vdw_lj_hbonded_hard);
  fprintf(outfile,"Use and extended vdW cutoff? (0:no,1:yes): %d\n",this.vdw_use_extended_cutoff);
  fprintf(outfile,"vdW function pointer allocated? (0:no,1:yes) %d",(this.vdw_function!=NULL));
  fprintf(outfile,"vdW clash energy at hard cutoff %g",this.vdw_clash_energy_at_hard_cutoff);
  fprintf(outfile,"extended vdW cutoff value %g\n",this.vdw_extended_cutoff);
  print_vdw_cutoff_distances(&this,outfile);
  /* stress */
  fprintf(outfile,"STRESS\n");
  fprintf(outfile,"stress_k %g\n",this.stress_k);
  fprintf(outfile,"stress_angle %g\n",this.stress_angle);

  /* hydrogen bond */
  fprintf(outfile,"H-BOND\n");
  fprintf(outfile,"OH cutoff square %g\n",this.hboh2);
  fprintf(outfile,"OHN angle cutoff %g\n",this.hbohn);
  fprintf(outfile,"COH angle cutoff %g\n",this.hbcoh);
  //fprintf(outfile,"OH cutoff decay width %g\n",this.hboh_decay_width);
  //fprintf(outfile,"OHN angle decay width %g\n",this.hbohn_decay_width);
  //fprintf(outfile,"COH angle decay width %g\n",this.hbcoh_decay_width);
  fprintf(outfile,"strength %g\n",this.hbs);
  /* contact parameters */
  fprintf(outfile,"CONTACT\n");
  fprintf(outfile,"cutoff square %g\n",this.touch2);
  fprintf(outfile,"part %g\n",this.part);
  fprintf(outfile,"split %g\n",this.split);
  fprintf(outfile,"sts %g\n",this.sts);
  /* biasing force constants */
  fprintf(outfile,"BIAS\n");
  fprintf(outfile,"contact map file %s\n",this.contact_map_file);
  fprintf(outfile,"eta_beta %g\n",this.bias_eta_beta);
  fprintf(outfile,"eta_alpha %g\n",this.bias_eta_alpha);
  fprintf(outfile,"kappa_alpha_3 %g\n",this.bias_kappa_alpha_3);
  fprintf(outfile,"kappa_alpha_4 %g\n",this.bias_kappa_alpha_4);
  fprintf(outfile,"kappa_beta %g\n",this.bias_kappa_beta);
  fprintf(outfile,"CA(0) or CB(1) %g\n",this.prt);
  fprintf(outfile,"r0_a %g\n",this.bias_r_alpha);
  fprintf(outfile,"r0_b %g\n",this.bias_r_beta);
  /* hydrophobicity */
  fprintf(outfile,"HYDROPHOBICITY\n");
  fprintf(outfile,"k_h %g\n",this.kauzmann_param);
  fprintf(outfile,"delta_decay %g\n",this.hydrophobic_cutoff_range);
  fprintf(outfile,"from i,i+%d\n",this.hydrophobic_min_separation);
  //fprintf(outfile,"hydrophobic_min_cutoff (if 1/dist) %g\n",this.hydrophobic_min_cutoff);
  //fprintf(outfile,"hydrophobic_max_cutoff (if 1/dist) %g\n",this.hydrophobic_max_cutoff);
  //fprintf(outfile,"hydrophobic_max_Eshift (if 1/dist) %g\n",this.hydrophobic_max_Eshift);
  //fprintf(outfile,"hydrophobic_r (if spline) %g\n",this.hydrophobic_r);
  //fprintf(outfile,"hydrophobic_half_delta (if spline) %g\n",this.hydrophobic_half_delta);
  /* electrostatics */
  fprintf(outfile,"ELECTROSTATICS\n");
  fprintf(outfile,"permittivity %g\n",this.recip_dielectric_param);
  fprintf(outfile,"Debye length %g\n",this.debye_length_param);
  fprintf(outfile,"from i,i+%d\n",this.electrostatic_min_separation);
  /* side chain hydrogen bond parameters */
  fprintf(outfile,"SIDE CHAIN H-BOND\n");
  fprintf(outfile,"H_s2b %g\n",this.sidechain_hbond_strength_s2b);
  fprintf(outfile,"H_b2s %g\n",this.sidechain_hbond_strength_b2s);
  fprintf(outfile,"H_s2s %g\n",this.sidechain_hbond_strength_s2s);
  fprintf(outfile,"cutoff %g\n",this.sidechain_hbond_cutoff);
  fprintf(outfile,"cos(angle cutoff) %g\n",this.sidechain_hbond_angle_cutoff);
  fprintf(outfile,"delta_decay %g\n",this.sidechain_hbond_decay_width);
  fprintf(outfile,"from i,i+%d\n",this.sidechain_hbond_min_separation);
  /* secondary radius of gyration */
  fprintf(outfile,"srgy_param %g\n",this.srgy_param);
  fprintf(outfile,"srgy_offset %g\n",this.srgy_offset);
  fprintf(outfile,"hphobic_srgy_param %g\n",this.hphobic_srgy_param);
  fprintf(outfile,"hphobic_srgy_offset %g\n",this.hphobic_srgy_offset);
  /* SS bond */
  fprintf(outfile,"S-bond strength %g\n", this.Sbond_strength);
  fprintf(outfile,"S-bond distance %g\n", this.Sbond_distance);
  fprintf(outfile,"S-bond cutoff %g\n",	this.Sbond_cutoff);
  fprintf(outfile,"S-bond dihedral cutoff %g\n",	this.Sbond_dihedral_cutoff);
  fprintf(outfile,"C-S-S angle hard coded\n");  
  fprintf(outfile,"FIXED AMINO ACIDS\n");
  fprintf(outfile,"fixed amino acid list file: %s\n",this.fixed_aalist_file);
  fprintf(outfile,"EXTERNAL CONSTRAINTS\n");
  fprintf(outfile,"constrained amino acid list file: %s\n",this.external_constrained_aalist_file);
  fprintf(outfile,"external type: %d\n",this.external_potential_type);
  fprintf(outfile,"external_k[0]: %g\n",this.external_k[0]);
  fprintf(outfile,"external_r0[0]: %g\n",this.external_r0[0]);
  if (this.external_potential_type == 3) fprintf(outfile,"external_ztip: %g\n",this.external_ztip);
  fprintf(outfile,"constrained amino acid list file(2): %s\n",this.external_constrained_aalist_file2);
  fprintf(outfile,"external type(2): %d\n",this.external_potential_type2);
  fprintf(outfile,"external_k[0](2): %g\n",this.external_k2[0]);
  fprintf(outfile,"external_r0[0](2): %g\n",this.external_r02[0]);
  if (this.external_potential_type2 == 3) fprintf(outfile,"external_ztip(2): %g\n",this.external_ztip2);
  fprintf(outfile,"============END=MODEL=PARAMETERS===================\n");

}


void flex_param_print(FLEX_params this, FILE *outfile){
  fprintf(outfile,"================FLEX=PARAMETERS===================\n");
  fprintf(outfile,"Number of processors: %d\n",this.number_of_processors);
  if(this.number_of_processors == 0) return;
  fprintf(outfile,"Modes: %d to %d\n",this.fromnode,this.tonode);
  fprintf(outfile,"Hbond strength for FLEX: %lf\n",this.hstrength_cutoff);
  if(this.only_bias_hbonds == 1.0){fprintf(outfile,"Using Only Bias H-bonds\n");}
  else{fprintf(outfile,"Using All H-bonds\n");}
  fprintf(outfile,"Max Totconf: %d Step: %lf\n",this.totalconf, this.step);
  fprintf(outfile,"Output path: %s FLEX dir: %s\n",this.output_path,this.flex_dir);
  fprintf(outfile, "Acceptance Rate Aim %lf\n",this.acceptance_rate_aim);
  fprintf(outfile,"FLEX command: %s\n",this.flex_cmd);
  fprintf(outfile,"============END=FLEX=PARAMETERS===================\n");
}

/* print simulation parameters (including model_parameters) */
void param_print(simulation_params this, FILE *outfile) {

  fprintf(outfile,"==============SIMULATION=PARAMETERS================\n");
  /* general simulation */
  fprintf(outfile,"--GENERAL--\n");
  fprintf(outfile,"infile_name %s\n",this.infile_name);
  fprintf(outfile,"outfile_name %s\n",this.outfile_name);
  fprintf(outfile,"pace %d\n",this.pace);
  fprintf(outfile,"stretch %d\n",this.stretch);
  fprintf(outfile,"test mask %x\n",this.tmask);
  fprintf(outfile,"random seed %d\n",this.seed);
  fprintf(outfile,"parameters %s\n",this.prm);
  fprintf(outfile,"acceptance ratio %g\n",this.acceptance_rate);
  fprintf(outfile,"amplitude %g\n",this.amplitude);
  fprintf(outfile,"keep_amplitude_fixed %d\n",this.keep_amplitude_fixed);
  fprintf(outfile,"acceptance %g\n",this.acceptance);
  fprintf(outfile,"accept_counter %d\n",this.accept_counter);
  fprintf(outfile,"reject_counter %d\n",this.reject_counter);
  fprintf(outfile,"acceptance_rate_tolerance %g\n",this.acceptance_rate_tolerance);
  fprintf(outfile,"amplitude_changing_factor %g\n",this.amplitude_changing_factor);
  if (this.MC_lookup_table != NULL && this.MC_lookup_table_n != NULL && this.Nchains > 0) {
    for (int i=0; i<4; i++) {
      for (int j=0; j<(this.NAA-1+this.Nchains);j++) {
        fprintf(outfile,"%d ",this.MC_lookup_table[i*(this.NAA+this.Nchains)+j]);
      }
      fprintf(outfile, " (%d valid elements)\n",this.MC_lookup_table_n[i]);
    }
  }
  /* peptide or multi-chain protein */
  fprintf(outfile,"--PEPTIDE--\n");
  fprintf(outfile,"number of chains: %d\n", this.Nchains);
  if (this.seq == NULL || this.sequence == NULL) {
    fprintf(outfile,"amino acid sequence %s\n",this.seq);
    fprintf(outfile,"amino acid sequence %s\n",this.sequence);
  } else {
    if (strlen(this.seq) == strlen(this.sequence)) { // 1 chain
      fprintf(outfile,"amino acid sequence %s\n",this.seq);
    } else {
      fprintf(outfile,"amino acid sequence %s\n",this.sequence);
      fprintf(outfile,"number of chains: %d\n", (int)strlen(this.sequence) - (int)strlen(this.seq));
    }
  }
  fprintf(outfile,"number of amino acids %d\n",this.NAA);
  /* thermodynamic beta */
  fprintf(outfile,"--TEMPERATURE--\n");
  fprintf(outfile,"thermobeta %g\n",this.thermobeta);
  fprintf(outfile,"lowtemp %d\n",this.lowtemp);
  fprintf(outfile,"beta1 %g\n",this.beta1);
  fprintf(outfile,"beta2 %g\n",this.beta2);
  fprintf(outfile,"bstp %g\n",this.bstp);
  fprintf(outfile,"intrvl %d\n",this.intrvl);
  /* nested sampling */
  fprintf(outfile,"--NESTED SAMPLING--\n");
  fprintf(outfile,"do nested sampling? %d\n",this.NS);
  fprintf(outfile,"NS iter %d\n",this.iter);
  fprintf(outfile,"NS iter_start %d\n",this.iter_start);
  fprintf(outfile,"NS iter_max %d\n",this.iter_max);
  fprintf(outfile,"NS logLstar %g\n",this.logLstar);
  fprintf(outfile,"NS logZ %g\n",this.logZ);
  fprintf(outfile,"NS logfactor %g\n",this.logfactor);
  fprintf(outfile,"NS alpha %g\n",this.alpha);
  fprintf(outfile,"NS logX %g\n",this.logX);
  fprintf(outfile,"NS logX_start %g\n",this.logX_start);
  fprintf(outfile,"NS Delta_logX %g\n",this.Delta_logX);
  fprintf(outfile,"NS log_DeltaX %g\n",this.log_DeltaX);
  fprintf(outfile,"NS H %g\n",this.H);
  fprintf(outfile,"NS N %d\n",this.N);
  /* checkpointing */
  fprintf(outfile,"--CHECKPOINTING--\n");
  fprintf(outfile,"num_NS_per_checkpoint %d\n",this.num_NS_per_checkpoint);
  fprintf(outfile,"checkpoint filename %s\n",this.checkpoint_filename);
  fprintf(outfile,"checkpoint_counter %d\n",this.checkpoint_counter);
  fprintf(outfile,"restart_from_checkpoint %d\n",this.restart_from_checkpoint);
  fprintf(outfile,"checkpoint %d\n",this.checkpoint);

  model_param_print(this.protein_model, outfile);
  flex_param_print(this.flex_params,outfile);
  fprintf(outfile,"==========END=SIMULATION=PARAMETERS================\n");

}
