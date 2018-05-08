/*
** Header with all the simulation parameters **
** and its default values
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/


#define DEFAULT_LONG_STRING_LENGTH 1024
#define DEFAULT_SHORT_STRING_LENGTH 256

#define NO_GAMMA 0 //no gamma atoms, only up to beta carbons
#define LINUS_GAMMA 1 //LINUS beta-gamma distances for all residues
#define CORRECT_GAMMA 2 //correct beta-gamma distances for all residues
#define CORRECT_KMQR_GAMMA 3 //LINUS beta-gamma distances for all but K,M,Q,R residues
#define DEFAULT_EXTENDED_VDW_CUTOFF_GAMMA 500
#define DEFAULT_EXTENDED_VDW_CUTOFF_NOGAMMA 100
#define HARD_CUTOFF_VDW_POTENTIAL 1
#define LJ_VDW_POTENTIAL 2

#define EXTERNAL_NEGATIVE 1001
#define EXTERNAL_NONE     1002
#define EXTERNAL_POSITIVE 1003
#define EXTERNAL_POSNEG   1004

//Usage help for the parameter string
#define PARAM_USE "Usage of the parameter string -p Param1=...[,...][,Param2=...][,Param3=...,...,...]\n\
Options:\n\
 Bias                   Bias potential\n\
 Hbond                  H-bond\n\
 Contact                contact \n\
 Gamma                  gamma atoms\n\
 VDW                    vdW potential\n\
 Stress                 backbone stress\n\
 Hydrophobic            hydrophobic\n\
 Fixit                  fix peptide\n\
 Elec                   electrostatic potential\n\
 SchHbond               side chain H-bond\n\
 Rgyr                   secondary radius of gyration\n\
 SSbond                 S-S bonds\n\
 fixed                  fixed amino acid list\n\
 external               external potential\n"

/* side chain properties of the protein model */
typedef struct {

  char id;                 // 1 letter code
  int number_of_gamma;     // number of gamma atoms in the amino acid
  /* geometry */
  double beta_gamma1_distance;                  // CB--G1 distance in Angstroms
  double alpha_beta_gamma1_angle;               // CA--CB--G1 angle in radians
  double beta_gamma2_distance;                  // CB--G2 distance in Angstroms
  double alpha_beta_gamma2_angle;               // CA--CB--G2 angle in radians
  double sidechain_dihedral_gauche_plus_prob;   // the probability of the N--CA--CB--G1 dihedral angle being gauche+
  double sidechain_dihedral_gauche_minus_prob;  // the probability of the N--CA--CB--G1 dihedral angle being gauche-
  /* vdW */
  double sidechain_vdw_radius1;        // exclusion radius of G1 atom in Angstroms
  double sidechain_vdw_radius2;        // exclusion radius of G2 atom in Angstroms
  double sidechain_vdw_depth1;         // exclusion depth of G1 atom in RT
  double sidechain_vdw_depth2;         // exclusion depth of G2 atom in RT
  double sidechain_vdw_depth1_sqrt;    // exclusion depth of G1 atom in RT sqrted
  double sidechain_vdw_depth2_sqrt;    // exclusion depth of G2 atom in RT sqrted
  /* hydrophobic interactions */
  int hydrophobic_atoms;                  // binary storage of which atoms are hydrophobic
  double hydrophobic_contact_radius_CB;   // hydrophobic contact radius of CB atom in Angstroms
  double hydrophobic_contact_radius_G1;   // hydrophobic contact radius of G1 atom in Angstroms
  double hydrophobic_contact_radius_G2;   // hydrophobic contact radius of G2 atom in Angstroms
  /* electrostatic interactions */
  int charged_atom;               // binary storage of which atom is charged (only 1 per side chain)
  double charged_atom_charge;     // charge on the charged atom
  /* H-bonds */
  int hydrogen_bond_donor_atoms;        // binary storage of which side chain atoms are hydrogen bond donors
  int hydrogen_bond_acceptor_atoms;     // binary storage of which side chain atoms are hydrogen bond acceptors
  double hydrogen_bond_donor_radius;    // size of side chain atoms as hydrogen bond donors
  double hydrogen_bond_acceptor_radius; // size of side chain atoms as hydrogen bond acceptors

} sidechain_properties_;

/* protein model and force field */
typedef struct {

  int fixit;              // fix peptide when reading it
  /* gamma atoms */
  int use_gamma_atoms;           // include gamma atoms in the model
  int use_original_gamma_atoms;  // use original sidechain angles when reading in PDB
  int use_3_states;              // adjust the sidechain angles to the stable states
				 // (+/-60,180 or +/-30 for PRO) when reading in PDB
  int fix_chi_angles;            // fix the read-in (or adjusted) side chain dihedral angles during the mc
  int fix_CA_atoms;              // fix the adjusted CA atomic positions (only rotate peptide bonds in mc)
  /* atomic radii */
  double rca;   // CA_ vdW Rmin
  double rcb;   // CB_ vdW Rmin
  double rc;    // C__ vdW Rmin
  double rn;    // N__ vdW Rmin
  double ro;    // O__ vdW Rmin
  double rs;    // S vdW Rmin
  double rring; // ring vdW Rmin
  double rel_vdw_cutoff;  // vdW cutoff relative to the Rmin
  int vdw_uniform_depth; // the same depth for all atoms
  double vdw_depth_ca; // CA_ vdW epsilon
  double vdw_depth_cb; // CB_ vdW epsilon
  double vdw_depth_c;  // C__ vdW epsilon
  double vdw_depth_n;  // N__ vdW epsilon
  double vdw_depth_o;  // O__ vdW epsilon
  double vdw_depth_s;  // S vdW epsilon
  double vdw_depth_ring; //ring vdW epsilon
  double vdw_depth_ca_sqrt; // CA_ vdW epsilon sqrt
  double vdw_depth_cb_sqrt; // CB_ vdW epsilon sqrt
  double vdw_depth_c_sqrt;  // C__ vdW epsilon sqrt
  double vdw_depth_n_sqrt;  // N__ vdW epsilon sqrt
  double vdw_depth_o_sqrt;  // O__ vdW epsilon sqrt
  double vdw_depth_s_sqrt;  // S vdW epsilon sqrt
  double vdw_depth_ring_sqrt;  // ring vdW epsilon sqrt
  // pairwise depths
  double vdw_depth_ca_ca;
  double vdw_depth_ca_cb;
  double vdw_depth_ca_c;
  double vdw_depth_ca_n;
  double vdw_depth_ca_o;
  double vdw_depth_cb_cb;
  double vdw_depth_cb_c;
  double vdw_depth_cb_n;
  double vdw_depth_cb_o;
  double vdw_depth_c_c;
  double vdw_depth_c_n;
  double vdw_depth_c_o;
  double vdw_depth_n_n;
  double vdw_depth_n_o;
  double vdw_depth_o_o;
  // pairwise energy shifts at the cutoff
  double vdw_shift;
  double vdw_Eshift_ca_ca;
  double vdw_Eshift_ca_cb;
  double vdw_Eshift_ca_c;
  double vdw_Eshift_ca_n;
  double vdw_Eshift_ca_o;
  double vdw_Eshift_cb_cb;
  double vdw_Eshift_cb_c;
  double vdw_Eshift_cb_n;
  double vdw_Eshift_cb_o;
  double vdw_Eshift_c_c;
  double vdw_Eshift_c_n;
  double vdw_Eshift_c_o;
  double vdw_Eshift_n_n;
  double vdw_Eshift_n_o;
  double vdw_Eshift_o_o;
  double vdw_backbone_cutoff; //the CA_ - CA_ distance cutoff for the backbone atoms up to CB
  double *vdw_gamma_gamma_cutoff; //the CA_ - CA_ distance cutoff for gamma - gamma vdW interactions 
  double *vdw_gamma_nongamma_cutoff; //the CA_ - CA_ distance cutoff for gamma - nongamma vdW interactions
  int vdw_use_extended_cutoff; //whether we should use an extended vdW cutoff, and its value
  double vdw_extended_cutoff; //the value of the extended vdW cutoff (default: 500 for models including gamma atoms, and 100 if gamma atoms are not used)
  int vdw_potential; //the potential form of the vdW interactions
  double (*vdw_function); // (vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift); //the vdW function
  double vdw_clash_energy_at_hard_cutoff; // energy penalty for clashing atoms
  int vdw_lj_neighbour_hard; //whether we should use the hard cutoff part of the LJ potential for those amino acids that are neighbours (useful for CD learning of vdW potential)
  int vdw_lj_hbonded_hard; // whether we should use the hard cutoff part of the LJ potential for those amino acids that are neighbours (useful for CD learning of vdW potential)

  /* stress */
  double stress_k; //force constant of stress
  double stress_angle; //the ideal N-CA-C angle

  /* hydrogen bond */
  double hboh2;  // cutoff distance square
  double hbohn;  // cos O-H-N angle cutoff
  double hbcoh;  // cos C-O-H angle cutoff
  //parameters for linear decay cutoff functions
  //double hboh_decay_width;  // cutoff distance decay width
  //double hbohn_decay_width;  // cos O-H-N angle cutoff decay width
  //double hbcoh_decay_width;  // cos C-O-H angle cutoff decay width
  double hbs;    // strength

  /* contact parameters */
  double touch2; // cutoff square
  double part;   // use CA(0) or CB(1)
  double split;  //
  double sts;    // strength

  /* biasing force constants */
  char *contact_map_file;    // contact map filename
  double bias_eta_beta;      // eta_beta
  double bias_eta_alpha;     // eta_alpha
  double bias_kappa_alpha_3; // kappa_alpha_3
  double bias_kappa_alpha_4; // kappa_alpha_4
  double bias_kappa_beta;    // kappa_beta
  double prt;                // use CA or CB
  double bias_r_alpha;       // r0_alpha
  double bias_r_beta;        // r0_beta

  /* hydrophobicity */
  double kauzmann_param;           // proportional to the Kauzmann parameter
  double hydrophobic_cutoff_range; // linear decay range for max->0 beyond the cutoff
  int hydrophobic_min_separation;  // minimum sequence separation
  // 1/dist potential form
  //double hydrophobic_min_cutoff;   // minimum distance cutoff for hydrophobic interactions, set to smaller than vdW cutoff
  //double hydrophobic_max_cutoff;   // maximum distance cutoff for hydrophobic interactions
  //double hydrophobic_max_Eshift;   // energy shift at the maximum distance cutoff (hydrophobic_max_cutoff^-1)
  // Spline potential form
  //double hydrophobic_r;            // the inflection point of the spline
  //double hydrophobic_half_delta;   // the distance over which the spline value changes by half the height of the potential

  /* electrostatics */
  double recip_dielectric_param;    // recip. permittivity
  double debye_length_param;        // Debye length
  int electrostatic_min_separation; // minimum sequence separation

  /* side chain hydrogen bond parameters */
  double sidechain_hbond_strength_s2b;  // sch -> bb h-bond strength
  double sidechain_hbond_strength_b2s;  // bb -> sch h-bond strength
  double sidechain_hbond_strength_s2s;  // sch -> sch h-bond strength
  double sidechain_hbond_cutoff;        // h-bond cutoff distance
  double sidechain_hbond_decay_width;   // linear decay range max->0 beyond the cutoff
  double sidechain_hbond_angle_cutoff;  // the cos(cutoff angle) for sidechain H-bonds
  int sidechain_hbond_min_separation;   // minimum sequence separation

  /* secondary radius of gyration */
  double srgy_param;
  double srgy_offset;
  double hphobic_srgy_param;
  double hphobic_srgy_offset;

  /* S-bond */
  double Sbond_strength;
  double Sbond_distance;
  double Sbond_cutoff;
  double Sbond_dihedral_cutoff;

  /* fixed amino acids */
  char *fixed_aalist_file;

  /* external potential */
  /* harmonic potential in the (+) and (-) directions with an offset */
  int external_potential_type;
  int external_direction[3]; /* - neither + both */
  double external_k[3]; /* x y z */
  double external_r0[3]; /* x y z */
  double external_ztip; /* for conincal potential */
  char *external_constrained_aalist_file;
  /* second potential */
  int external_potential_type2;
  int external_direction2[3]; /* - neither + both */
  double external_k2[3]; /* x y z */
  double external_r02[3]; /* x y z */
  double external_ztip2; /* for conincal potential */
  char *external_constrained_aalist_file2;
  /* sidechain properties */
  sidechain_properties_ *sidechain_properties;

} model_params;


typedef struct {
  int number_of_processors;
  char *output_path;
  char *outputpdb_filename;
  char *flex_dir;
  char **filenames_to_read_in;
  int size_of_filename_to_read_in;
  double hstrength_cutoff;
  int fromnode, tonode, freq, totalconf;
  int only_bias_hbonds; //= 0 if not, 1 if
  char *flex_cmd;
  double step;
  double acceptance_rate_aim;
  double acceptance_rate_tolerance;
  double flex_changing_factor;
  int MCiter;
  int mode;
} FLEX_params ;

/* simulation parameters */
typedef struct {

  /* general simulation */
  char *infile_name;    /* input filename */
  char *outfile_name;   /* output filename */
  FILE *infile;
  FILE *outfile;
  unsigned int pace;
  unsigned int stretch;
  unsigned int tmask;
  unsigned int seed;
  char *prm;
  double acceptance_rate;
  double amplitude;
  int keep_amplitude_fixed;
  double acceptance;
  int accept_counter;
  int reject_counter;
  double acceptance_rate_tolerance;
  double amplitude_changing_factor;
  int number_initial_MC; //number of initial MC moves before NS starts
  int *MC_lookup_table; //lookup table for random moves
  int *MC_lookup_table_n; //the number of valid elements for each loop length

//  char *infile; //use stdin
  /* peptide */
  char *seq;
  char *sequence;	//sequence with '_' separators for chain breaks
  int NAA;		// not used
  int Nchains;
  /* thermodynamic beta (reciprocal kT) schedule for annealing or tempering */
  double thermobeta;
  int lowtemp;
  double beta1;
  double beta2;
  double bstp;
  unsigned int intrvl;

  /* various tests */
  double *energy_gradient;
  double *energy_probe_1_this;
  double *energy_probe_1_last;
  int *energy_probe_1_calc;

  /* nested sampling */
  int NS;          /* whether NS calculation is turned on(=1) or not (0) */
  int iter;        /* iteration counter */
  int iter_start;  /* starting value of iteration */
  int iter_max;    /* end value of iteration */
  double logLstar; /* log L* */
  double logZ;     /* log evidence */
  double logfactor;/* log factor for parallel nested sampling */
  double H;        /* information */
  int N;           /* the number of NS points */

  double alpha;    //shrinkage ratio of X, the available prior space ratio
  double logX;     //log of the available prior space ratio estimate
  double logX_start; //start value of log of the available prior space ratio estimate (important for restarts!)
  double Delta_logX; //step length along the logX axis: logX_{i-1} - logX_i = log alpha
  double log_DeltaX; //current width of the prior space weight: log( X_{i-1} - X_i ) = i * log(1 - alpha)

  /* checkpointing */
  int num_NS_per_checkpoint;
  char *checkpoint_filename;
  FILE *checkpoint_file;
  int checkpoint_counter;
  int checkpoint;
  int restart_from_checkpoint;
  model_params protein_model;

  /*normal mode analysis parameters */
  FLEX_params flex_params;

} simulation_params;


/* initialising and finalising */
void param_initialise(simulation_params *this);
void param_finalise(simulation_params *this);
void model_param_initialise(model_params *this);
void model_param_finalise(model_params *this);

void flex_params_initialise(FLEX_params *this);
void flex_params_finalise(FLEX_params *this);
void flex_setup_command(FLEX_params *this);

void vdw_param_zero(model_params *this);

/* CD learnt default parameters */
void set_lj_default_params(model_params *this);
void set_hard_cutoff_default_params(model_params *this);

/* copying */
void sim_params_copy(simulation_params *to, simulation_params *from);
void model_params_copy(model_params *to, model_params *from);
void flex_params_copy(FLEX_params *to, FLEX_params *from);
void copy_string(char** const to, const char* const from);

/* reading in and processing input */
void model_param_read(char *prm, model_params *this,FLEX_params *nma_params);
int flex_param_read(char *prm, FLEX_params *this);
void vdw_param_calculate(model_params *this);

/* printing */
void print_vdw_cutoff_distances(model_params *mod_params, FILE *outfile);
void model_param_print(model_params this, FILE *outfile);
void param_print(simulation_params this, FILE *outfile);
void flex_param_print(FLEX_params this, FILE *outfile);
