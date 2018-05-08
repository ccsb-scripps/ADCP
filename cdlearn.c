/*
** ML inference of model parameters using Contrastive Divergence.
**
** Copyright (c) 2004 - 2008 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include<signal.h>
#include<math.h>
#include<omp.h>

#include"error.h"
#include"params.h"
#include"aadict.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"
#include"vdw.h"
#include"energy.h"
#include"metropolis.h"
#include"probe.h"
#include"flex.h"
#include"cdlearn.h"

#define VER "CD learn program of CRANKITE (c)\n"
#define USE "Usage: %s [options]\n\
Options:\n\
 -o outfile           output file of CD iterations\n\
 -l learn_string      the parameters to be learnt(VWTHBRKDESGF)\n\
 -L list_file         the list file with the names of PDBs and ICMs without the extension\n\
 -i max_iter          The number of iterations for the CD learning\n\
 -K n_mc_steps        number of MC steps for generating the perturbed distribution\n\
 -A AMPLITUDE,FIX_AMP crankshaft rotation amplitude and whether the amplitude should be kept fixed (default is no fixing) \n\
 -R name,N            restart file name to read initial parameters from, the restart files will be filename_* at every 100 iterations\n\
 -I N                 restart iteration count, the restart file will be restartfilename_N\n\
 -p STRING            undocumented custom parameters\n\
 -d vdw_model         Which vdW model to set the defaults for (default: LJ, could also be hard_cutoff)\n\
 -a DOUBLE            if a>0, use adaptive learning rate\n\
 -s SEED              random seed\n\
 -t MASK,OPTIONS      hexadecimal mask of active tests\n"


/* K (=sim_params->pace) steps long MCMC for Contrastive Divergence learning.
   Possibly recalculate the amplitude to result the required acceptance rate,
   but do it less and less frequently, as it converges. */
double simulate(Chain * chain, Chaint *chaint, Biasmap* biasmap, simulation_params *sim_params, int i_protein, int iter)
{
	double temp; //unused for MC
	
	Chain *chain2 = (Chain *)malloc(sizeof(Chain));
	chain2->aa = NULL; chain2->xaa = NULL; chain2->erg = NULL; chain2->xaa_prev = NULL;
	allocmem_chain(chain2,chain->NAA,chain->Nchains);
	
	/* do tests at step 0 */
	tests(chain,biasmap,sim_params->tmask, sim_params, 0x11, NULL);
	double old_energy = totenergy(chain);

	/* possibly adjust the amplitude for the acceptance rate */
	if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
	    if (iter<=50 || (iter > 50 && iter <= 200 && iter % 10 == 0) || (iter > 200 && iter % 25 == 0)) {
		if (iter > 50) sim_params->amplitude_changing_factor = 0.95;
		if (iter > 200) sim_params->amplitude_changing_factor = 0.98;
		//fprintf(stderr,"Iter %d acceptance tolerance %g amplitude changing factor %g\n",iter,sim_params->acceptance_rate_tolerance,sim_params->amplitude_changing_factor);
		//fprintf(stderr,"Acceptance rate %d %g\n",i_protein,sim_params->acceptance_rate);
		/* This bit ensures the amplitude of the moves
		is independent of the chain's history*/
		copybetween(chain2,chain);
		move(chain2,chaint,biasmap,0.0,&temp,-1, sim_params);
		for(int j = 1; (j < sim_params->pace || j < 1024); j++){
		    move(chain2,chaint,biasmap,0.0,&temp,1, sim_params);	
		}
	    }
	    fprintf(stderr,"Amplitude %d %g\n",i_protein,sim_params->amplitude);
	}

	/* take the K steps of MC for the contrastive divergence */
	for (int j = 0; j < sim_params->pace; j++) {
	    move(chain,chaint,biasmap,0.0,&temp,0, sim_params);
	}

	/* do tests at K steps */
	tests(chain,biasmap,sim_params->tmask, sim_params, 0x11, NULL);

	double new_energy = totenergy(chain);
	//fprintf(stderr,"Energy_change %d %g %g %g\n",i_protein,old_energy,new_energy,new_energy-old_energy);
#ifdef DEBUG
	char my_string[DEFAULT_LONG_STRING_LENGTH] = "";
	sprintf(my_string,"%d ",i_protein);
	for (int i=0; i<34; i++) {
	    sprintf(my_string,"%s %f", my_string, sim_params->energy_probe_1_this[i] - sim_params->energy_probe_1_last[i]);
	}
	sprintf(my_string,"%s\n", my_string);
#endif

	freemem_chain(chain2); free(chain2);

	return new_energy - old_energy;
	
}

/* Read in all sorts of parameters for the CD learning. */
char *read_options(int argc, char *argv[], simulation_params *sim_params, cdlearn_params *cd_params)
{
	unsigned int pace = 0, tmask = 0x0, seed = 0;
	int iter_max = 1000;

	int i, opt;
	char *retval = NULL;
	char *learn_string = NULL;
	char *restart_filename = NULL;
	char *list_name = NULL;
	double adaptive = cd_params->adaptive_learning_param;
	double amplitude = sim_params->amplitude;
	int keep_amplitude_fixed = sim_params->keep_amplitude_fixed;
	char error_string[DEFAULT_LONG_STRING_LENGTH]="";


	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			if (freopen(argv[i], "r", stdin) == NULL)
				retval = argv[i];
			continue;
		}

		opt = argv[i][1];
		if (++i >= argc && opt != 'n')
			opt = 0;

		switch (opt) {
		case 'i':
			sscanf(argv[i], "%u", &iter_max);
			cd_params->iter_max = iter_max;
			break;
		case 'o':
			copy_string(&(cd_params->output_filename),argv[i]);
			break;
		case 'p':
			copy_string(&(sim_params->prm),argv[i]);
			break;
		case 'd':
			if (strcmp(argv[i],"hard_cutoff")==0) {
				sim_params->protein_model.vdw_potential=HARD_CUTOFF_VDW_POTENTIAL;
				set_hard_cutoff_default_params(&(sim_params->protein_model));
			} else if (strcmp(argv[i],"lj")==0) {
				sim_params->protein_model.vdw_potential=LJ_VDW_POTENTIAL;
				set_lj_default_params(&(sim_params->protein_model));
			} else {
				sprintf(error_string,"Unknown value for vdW potential model (%s).  It must be one of hard_cutoff and lj.",argv[i]);
				stop(error_string);
			}
			break;
		case 'K':
			sscanf(argv[i], "%u", &pace);
			sim_params->pace = pace;
			sim_params->stretch = 1;
			break;
		case 'R':
			copy_string(&restart_filename,argv[i]);
			cd_params->restart_filename = realloc(cd_params->restart_filename,DEFAULT_LONG_STRING_LENGTH);
			strcpy(cd_params->restart_filename,restart_filename);
			fprintf(stderr,"The restart filename is %s",cd_params->restart_filename);
			free(restart_filename);
			break;
		case 'I':
			sscanf(argv[i], "%d", &cd_params->iter_start);
			fprintf(stderr,"The restart iter count is %d",cd_params->iter_start);
			break;
		case 's':
			sscanf(argv[i], "%u", &seed);
			sim_params->seed = seed;
			break;
		case 'l':
			copy_string(&learn_string,argv[i]);
			cd_params->learn_string = realloc(cd_params->learn_string,DEFAULT_LONG_STRING_LENGTH);
			strcpy(cd_params->learn_string,learn_string);
			free(learn_string);
			fprintf(stderr,"The learn string is %s",cd_params->learn_string);
			break;
		case 'L':
			copy_string(&list_name,argv[i]);
			cd_params->pdb_list_filename = realloc(cd_params->pdb_list_filename,DEFAULT_LONG_STRING_LENGTH);
			strcpy(cd_params->pdb_list_filename,list_name);
			free(list_name);
			break;
		case 't':
			sscanf(argv[i], "%x", &tmask);
			sim_params->tmask = tmask;
			break;
		case 'a':
			sscanf(argv[i], "%lf", &adaptive);
			cd_params->adaptive_learning_param = adaptive;
			break;
		case 'A': //amplitude
			sscanf(argv[i], "%lf,%d", &amplitude,&keep_amplitude_fixed);
			if(amplitude > 0.0 || amplitude <= -M_PI){
			  fprintf(stderr,"Amplitude must be between -PI and 0.0, setting it to -0.01\n");
			  amplitude = -0.01;
			}
			sim_params->amplitude = amplitude;
			sim_params->keep_amplitude_fixed = keep_amplitude_fixed;
			break;
		default:
			fprintf(stderr, VER USE PARAM_USE, argv[0]);
			helps();
			exit(EXIT_FAILURE);
		}
	}

	if (retval == NULL) {
	  sim_params->sequence=NULL;
	  sim_params->seq=NULL;
	} else {
	  sim_params->sequence = realloc(sim_params->sequence,strlen(retval+1));
	  strcpy(sim_params->sequence,retval);
	  sim_params->seq = realloc(sim_params->seq,strlen(retval+1));
	  strcpy(sim_params->seq,retval);
	}
	return retval;
}

void graceful_exit(int sig)
{
	exit(EXIT_FAILURE);	/* flushes streams too */
}


void set_random_seed(simulation_params *sim_params) {

	if (sim_params->seed == 0)
		sim_params->seed = (unsigned int) time(NULL);

	/*random seed*/
	srand(sim_params->seed);

}


void update_sim_params_from_cd_learn(simulation_params *sim_params, cdlearn_params *cd_params) {

  /* WHAT GRADIENT TO CALCULATE */
  /* only those that would change */
  int i;
  for (i=0; i<36; i++) {
    if (cd_params->raa[i] == 0) {
	sim_params->energy_probe_1_calc[i] = 0;
    } else {
	sim_params->energy_probe_1_calc[i] = 1;
    }
  }

  /* PARAMETERS */
  /* hydrogen bond */
  sim_params->protein_model.hboh2 = cd_params->aa[0];
  sim_params->protein_model.hbohn = cd_params->aa[1];
  sim_params->protein_model.hbcoh = cd_params->aa[2];
  sim_params->protein_model.hbs = cd_params->aa[3];
  /* biasing force constants */
  sim_params->protein_model.bias_eta_beta = cd_params->aa[4];
  sim_params->protein_model.bias_eta_alpha = cd_params->aa[5];
  sim_params->protein_model.bias_kappa_alpha_3 = cd_params->aa[6];
  sim_params->protein_model.bias_kappa_alpha_4 = cd_params->aa[7];
  sim_params->protein_model.bias_kappa_beta = cd_params->aa[8];
  sim_params->protein_model.bias_r_alpha = cd_params->aa[12];
  sim_params->protein_model.bias_r_beta = cd_params->aa[13];
  /* hydrophobicity */
  sim_params->protein_model.kauzmann_param = cd_params->aa[9];
  sim_params->protein_model.hydrophobic_cutoff_range = cd_params->aa[10];
  //sim_params->protein_model.hydrophobic_min_separation = 2;
  /* electrostatics */
  sim_params->protein_model.recip_dielectric_param = cd_params->aa[11];
  sim_params->protein_model.debye_length_param = cd_params->aa[14];
  //sim_params->protein_model.electrostatic_min_separation = cd_params->aa[0];
  /* side chain hydrogen bond parameters */
  sim_params->protein_model.sidechain_hbond_strength_s2b = cd_params->aa[15];
  sim_params->protein_model.sidechain_hbond_strength_b2s = cd_params->aa[16];
  sim_params->protein_model.sidechain_hbond_strength_s2s = cd_params->aa[17];
  //sim_params->protein_model.sidechain_hbond_cutoff = cd_params->aa[0];
  //sim_params->protein_model.sidechain_hbond_decay_width = cd_params->aa[0];
  sim_params->protein_model.sidechain_hbond_angle_cutoff = cd_params->aa[18];
  //sim_params->protein_model.sidechain_hbond_min_separation = cd_params->aa[0];
  /* atomic radii */
  sim_params->protein_model.rca = cd_params->aa[21];
  sim_params->protein_model.rcb = cd_params->aa[22];
  sim_params->protein_model.rc = cd_params->aa[23];
  sim_params->protein_model.rn = cd_params->aa[24];
  sim_params->protein_model.ro = cd_params->aa[25];
  sim_params->protein_model.rs = cd_params->aa[26];
  sim_params->protein_model.vdw_depth_ca = cd_params->aa[28];
  sim_params->protein_model.vdw_depth_cb = cd_params->aa[29];
  sim_params->protein_model.vdw_depth_c = cd_params->aa[30];
  sim_params->protein_model.vdw_depth_n = cd_params->aa[31];
  sim_params->protein_model.vdw_depth_o = cd_params->aa[32];
  sim_params->protein_model.vdw_depth_s = cd_params->aa[33];
  /* stress */
  sim_params->protein_model.stress_k = cd_params->aa[27];
  /* contact parameters */
  //sim_params->protein_model.touch2 = 49;
  //sim_params->protein_model.part = cd_params->aa[0];
  //sim_params->protein_model.split = cd_params->aa[0];
  //sim_params->protein_model.sts = cd_params->aa[0];
  /* secondary radius of gyration */
  sim_params->protein_model.srgy_param = cd_params->aa[19];
  sim_params->protein_model.srgy_offset = cd_params->aa[20];
  sim_params->protein_model.hphobic_srgy_param = cd_params->aa[34];
  sim_params->protein_model.hphobic_srgy_offset = cd_params->aa[35];
  /* S-S bond */
  //sim_params->protein_model.Sbond_strength = cd_params->aa[0];
  //sim_params->protein_model.Sbond_distance = cd_params->aa[0];
  //sim_params->protein_model.Sbond_cutoff = cd_params->aa[0];
  //sim_params->protein_model.Sbond_dihedral_cutoff = cd_params->aa[0];

  vdw_param_calculate(&(sim_params->protein_model));

}

/* Initialise the CD learning parameters to default values (mainly 0s)
   and then reading in initial values of CD learning parameters from a restart file */
void cd_learn_param_initialise(cdlearn_params *cd_params){

  /* default */
  if (cd_params->iter_max == 0) cd_params->iter_max = 2;
  int i;
  for (i=0; i<36; i++) {
     cd_params->aa[i] = 0;
     cd_params->raa[i] = 0;
  }
  cd_params->aa[21]=1.85;
  cd_params->aa[22]=2.00;
  cd_params->aa[23]=1.85;
  cd_params->aa[24]=1.75;
  cd_params->aa[25]=1.60;
  cd_params->aa[26]=2.00;
  cd_params->aa[28]=0.2;
  cd_params->aa[29]=0.2;
  cd_params->aa[30]=0.2;
  cd_params->aa[31]=0.2;
  cd_params->aa[32]=0.2;
  cd_params->aa[33]=0.2;
  cd_params->aa[27]=70.0;
  cd_params->aa[10]=2.8;


  /* default from string */
  fprintf(stderr,"Setting up default from string: %s\n",cd_params->learn_string);
  if (cd_params->learn_string==NULL) stop("No learn string was given.  Abort.");
  char learn_what;
  for (i=0; i< strlen(cd_params->learn_string); i++) {
     learn_what = cd_params->learn_string[i];
     fprintf(stderr,"%c",learn_what);
     switch(learn_what) {
	case 'H': cd_params->aa[0]=4.25; cd_params->raa[0]=.002;
		  cd_params->aa[1]=0.89; cd_params->raa[1]=.0002;
		  cd_params->aa[2]=0.77; cd_params->raa[2]=.0002;
		  cd_params->aa[3]=4.60; cd_params->raa[3]=0.01;
		  break;
	case 'B': cd_params->aa[4]=3.8; cd_params->raa[4]=.5;
		  cd_params->aa[5]=10.6; cd_params->raa[5]=1.;
		  //cd_params->aa[6]=2.5; cd_params->raa[6]=.005
		  cd_params->aa[8]=1.3; cd_params->raa[8]=.005;
		  break;
	case 'K': cd_params->aa[9]=0.08; cd_params->raa[9]=0.0001;
		  cd_params->aa[10]=2.8; cd_params->raa[10]=0;
		  break;
	case 'D': cd_params->aa[11]=11.5; cd_params->raa[11]=1.0;
		  break;
	case 'E': cd_params->aa[14]=2.5; cd_params->raa[14]=0.1;
		  break;
	case 'R': cd_params->aa[12]=5.45; cd_params->raa[12]=0.01;
		  cd_params->aa[13]=5.4; cd_params->raa[13]=0.01;
		  break;
	case 'S': cd_params->aa[15]=1.0; cd_params->raa[15]=0.1;
		  cd_params->aa[16]=1.0; cd_params->raa[16]=0.1;
		  cd_params->aa[17]=1.0; cd_params->raa[17]=0.1;
		  cd_params->aa[18]=-0.5; cd_params->raa[18]=0.01;
		  break;
	case 'G': cd_params->aa[19]=1.; cd_params->raa[19]=0.1;
		  cd_params->aa[20]=5.; cd_params->raa[20]=10.;
		  break;
	case 'F': cd_params->aa[34]=1.; cd_params->raa[34]=0.1;
		  cd_params->aa[35]=3.; cd_params->raa[35]=10.;
		  break;
	case 'V': cd_params->aa[21]=2.00; cd_params->raa[21]=0.0020;
		  cd_params->aa[22]=2.00; cd_params->raa[22]=0.0001;
		  cd_params->aa[23]=1.85; cd_params->raa[23]=0.0001;
		  cd_params->aa[24]=2.00; cd_params->raa[24]=0.0001;
		  cd_params->aa[25]=1.60; cd_params->raa[25]=0.0001;
		  cd_params->aa[26]=2.00; cd_params->raa[26]=0.0010;
		  break;
	case 'W': cd_params->aa[28]=0.2; cd_params->raa[28]=0.0001;
		  cd_params->aa[29]=0.2; cd_params->raa[29]=0.0001;
		  cd_params->aa[30]=0.2; cd_params->raa[30]=0.0002;
		  cd_params->aa[31]=0.2; cd_params->raa[31]=0.0002;
		  cd_params->aa[32]=0.2; cd_params->raa[32]=0.0002;
		  cd_params->aa[33]=0.2; cd_params->raa[33]=0.0010;
		  break;
	case 'T': cd_params->aa[27]=70.; cd_params->raa[27]=10.;
		  break;
     }
  }
  fprintf(stderr,"\n");

  /* initialise from restart file if it exists */
  if (cd_params->restart_filename == NULL) {
     fprintf(stderr,"WARNING! No restart file given.\n");
     if (strlen(cd_params->learn_string)==0) stop("The learn string is empty and no restart file given.");
  } else {
     FILE *restart_file;
     char *restart_filename = NULL;
     restart_filename = (char *)realloc(restart_filename,(strlen(cd_params->restart_filename)+10)*sizeof(char));
     if (cd_params->iter_start==0) {
         strcpy(restart_filename,cd_params->restart_filename);
     } else {
         sprintf(restart_filename,"%s_%d",cd_params->restart_filename,cd_params->iter_start);
     }
     restart_file = fopen(restart_filename,"r");
     int k;

     if (restart_file == NULL) {
        fprintf(stderr,"WARNING! No restart file given.\n");
        if (strlen(cd_params->learn_string)==0) stop("The learn string is empty and no restart file given.");
     } else {
        fprintf(stderr,"Setting CD learn initial parameters from restart file %s.",restart_filename);
        k = fscanf(restart_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&cd_params->aa[0], &cd_params->aa[1], &cd_params->aa[2], &cd_params->aa[3], &cd_params->aa[4],
		&cd_params->aa[5], &cd_params->aa[6], &cd_params->aa[7], &cd_params->aa[8], &cd_params->aa[9],
		&cd_params->aa[10], &cd_params->aa[11], &cd_params->aa[12], &cd_params->aa[13], &cd_params->aa[14],
		&cd_params->aa[15], &cd_params->aa[16], &cd_params->aa[17], &cd_params->aa[18], &cd_params->aa[19],
		&cd_params->aa[20], &cd_params->aa[21], &cd_params->aa[22], &cd_params->aa[23], &cd_params->aa[24],
		&cd_params->aa[25], &cd_params->aa[26], &cd_params->aa[27], &cd_params->aa[28], &cd_params->aa[29],
		&cd_params->aa[30], &cd_params->aa[31], &cd_params->aa[32], &cd_params->aa[33], &cd_params->aa[34],
		&cd_params->aa[35]);

        if (k != 36) {
	   fclose(restart_file);
	   restart_file = NULL;
	   stop("Problem while reading aa-s from restart file.");
        }

        k = fscanf(restart_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&cd_params->raa[0], &cd_params->raa[1], &cd_params->raa[2], &cd_params->raa[3], &cd_params->raa[4],
		&cd_params->raa[5], &cd_params->raa[6], &cd_params->raa[7], &cd_params->raa[8], &cd_params->raa[9],
		&cd_params->raa[10], &cd_params->raa[11], &cd_params->raa[12], &cd_params->raa[13], &cd_params->raa[14],
		&cd_params->raa[15], &cd_params->raa[16], &cd_params->raa[17], &cd_params->raa[18], &cd_params->raa[19],
		&cd_params->raa[20], &cd_params->raa[21], &cd_params->raa[22], &cd_params->raa[23], &cd_params->raa[24],
		&cd_params->raa[25], &cd_params->raa[26], &cd_params->raa[27], &cd_params->raa[28], &cd_params->raa[29],
		&cd_params->raa[30], &cd_params->raa[31], &cd_params->raa[32], &cd_params->raa[33], &cd_params->raa[34],
		&cd_params->raa[35]);

        if (k != 36) {
	   fclose(restart_file);
	   restart_file = NULL;
	   stop("Problem while reading raa-s from restart file.");
        }
        fclose(restart_file);
	restart_file = NULL;
     }
     free(restart_filename);
  }

}

/* Finalise the CD learning parameters */
void cd_param_finalise(cdlearn_params *cd_params){
	if (cd_params->learn_string != NULL) free(cd_params->learn_string);
	if (cd_params->restart_filename != NULL) free(cd_params->restart_filename);
	if (cd_params->pdb_list_filename != NULL) free(cd_params->pdb_list_filename);
	if (cd_params->output_filename != NULL) free(cd_params->output_filename);
}

/* Write restart file with CD learning parameters */
void cd_learn_write_restart_file(cdlearn_params *this, char *restart_filename){

  FILE *restart_file;

  restart_file = fopen(restart_filename,"w");

  fprintf(restart_file,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
		this->aa[0], this->aa[1], this->aa[2], this->aa[3], this->aa[4],
		this->aa[5], this->aa[6], this->aa[7], this->aa[8], this->aa[9],
		this->aa[10], this->aa[11], this->aa[12], this->aa[13], this->aa[14],
		this->aa[15], this->aa[16], this->aa[17], this->aa[18], this->aa[19],
		this->aa[20], this->aa[21], this->aa[22], this->aa[23], this->aa[24],
		this->aa[25], this->aa[26], this->aa[27], this->aa[28], this->aa[29],
		this->aa[30], this->aa[31], this->aa[32], this->aa[33], this->aa[34],
		this->aa[35]);

  fprintf(restart_file,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
		this->raa[0], this->raa[1], this->raa[2], this->raa[3], this->raa[4],
		this->raa[5], this->raa[6], this->raa[7], this->raa[8], this->raa[9],
		this->raa[10], this->raa[11], this->raa[12], this->raa[13], this->raa[14],
		this->raa[15], this->raa[16], this->raa[17], this->raa[18], this->raa[19],
		this->raa[20], this->raa[21], this->raa[22], this->raa[23], this->raa[24],
		this->raa[25], this->raa[26], this->raa[27], this->raa[28], this->raa[29],
		this->raa[30], this->raa[31], this->raa[32], this->raa[33], this->raa[34],
		this->raa[35]);

  fclose(restart_file);
  restart_file = NULL;

}


/* If the values are meaning less, e.g. VDW radii are negative, adjust them */
void correct_negative_aa_values(cdlearn_params *cd_params) {

    for (int k=0; k<36; k++) {
	if (k==18) continue;
	if (cd_params->aa[k] < 0) {
	    fprintf(stderr,"WARNING!  cd_params->aa[%d] had to be adjusted!  It was %g, now 0.\n",k,cd_params->aa[k]);
	    cd_params->aa[k] = 0;
	}
    }

}


/* contact map printing using OXZU characters */
void print_contact_map(double *distb, int NAA, FILE *outfile) {
			
	int i, j;
	char c;

	for (i = 1; i < NAA; i++) {
		for (j = 1; j < NAA; j++) {
			int num = nearbyint(distb[i * (NAA ) + j]);
			switch (num) {
			case 0:
				c=*"O";
				break;
			case -1:
				c=*"Z";
				break;
			case 1:
				if (abs(i-j)<=1) {
					c=*"U";
				} else {
					c=*"X";
				}
				break;
			default:
				fprintf(stderr,"Invalid contact map value %g\n",distb[i * (NAA + 1) + j]);
				exit(EXIT_FAILURE);
			}
			fputc(c,outfile);
//			fprintf(outfile,"%g ",distb[i * NAA + j]);
		}
		fputc('\n',outfile);
	}

}


int main(int argc, char *argv[])
{

	time_t timer1 = time(NULL);

	/* print argument line */
	for (int i=0; i<argc; i++) {
		fprintf(stderr,"%s ",argv[i]);
	}
	fprintf(stderr,"\n");

	char *seq;
	simulation_params sim_params;
	cdlearn_params cd_params;
	cd_params.learn_string = NULL;
	cd_params.restart_filename = NULL;
	cd_params.iter_start = 0;
	cd_params.pdb_list_filename = NULL;
	cd_params.output_filename = NULL;
	cd_params.iter_max = 0;
	cd_params.adaptive_learning_param = 0;
	for (int i=0; i<36; i++) {
	    cd_params.total_energy_gradient[i] = 0;
	}

	signal(SIGTERM, graceful_exit);


	/* SET THE MAIN SIMULATION PARAMS */
	param_initialise(&sim_params); //set default
	sim_params.protein_model.vdw_uniform_depth = 1; // use the same vdW depth for all atoms by default!
	fprintf(stderr,"Using a uniform vdW depth by default (may be overwritten in -p VDW=...).");
	sim_params.protein_model.use_gamma_atoms = CORRECT_KMQR_GAMMA; // use the correct beta-gamma distances for the K, M, Q and R residues, and the LINUS distances for all other residues
	fprintf(stderr,"Using LINUS distances for all but K, M, Q and R residues by default (may be overwritten in -p Gamma=...).");
	sim_params.protein_model.vdw_use_extended_cutoff = 1; // use the extended vdW cutoff distance
	fprintf(stderr,"Using the extended cutoff (may be overwritten in -p VDW=...).");
	sim_params.protein_model.vdw_potential = LJ_VDW_POTENTIAL;
	fprintf(stderr,"Using the LJ potential as default (may be overwritten in -p VDW=...).");
	/* default crankshaft max amplitude as 0.01 radian, fixed */
	sim_params.amplitude = -0.01;
	sim_params.keep_amplitude_fixed = 1;
	//param_print(sim_params, stdout); //default

	seq = read_options(argc, argv, &sim_params, &cd_params);

	/* setting default parameters */
	if(sim_params.tmask==0x0) {
	    sim_params.tmask = 0x2000;
	    fprintf(stderr,"Setting default test mask %x",sim_params.tmask);
	} else {
	    sim_params.tmask |= 0x2000;
	}
	if(sim_params.pace==0) {
	    fprintf(stderr,"Setting default pace and stretch 4096x1");
	    sim_params.pace = 4096; sim_params.stretch = 1;
	}
	if(cd_params.iter_max == 0) {
	    fprintf(stderr,"Setting iteration limit to default: 1000");
	    cd_params.iter_max = 1000;
	}
	/* temperature of simulation (room temperature) */
	sim_params.thermobeta = 1.0;
	/* random seed */
	set_random_seed(&sim_params);


	/* INITIALISE PROTEIN MODEL WITHIN MAIN SIMULATION PARAMETERS */
	model_param_read(sim_params.prm,&(sim_params.protein_model),&(sim_params.flex_params));
	initialize_sidechain_properties(&(sim_params.protein_model));
	vdw_cutoff_distances_calculate(&sim_params, stderr, 0);
	peptide_init();
	//param_print(sim_params, stdout); //read-in


	/* READ IN ALL PDB CHAINS */
	fprintf(stderr,"Creating PDB library.\n");
	/* allocate memory for the original PDB library */
	Chain *temporary = (Chain *)malloc(sizeof(Chain));
	temporary->NAA = 0; temporary->aa = NULL; temporary->xaa = NULL; temporary->erg = NULL; temporary->xaa_prev = NULL;

	/* PDB library */
	Chain *all_chains = NULL;
	Chain *all_chains_sim = NULL;
	Chaint *all_chaints = NULL;
	Biasmap *all_biasmaps = NULL;
 
	FILE *list_file;
	if (cd_params.pdb_list_filename == NULL) stop("No list file given (-L option).");
	list_file = fopen(cd_params.pdb_list_filename,"r");
	if (list_file == NULL) stop("List file does not exist (-L option).");

	int n_proteins = 0;
	char next_pdb_filename[DEFAULT_LONG_STRING_LENGTH];
	char pdb_filename[DEFAULT_LONG_STRING_LENGTH];
	while (fscanf(list_file, "%s", next_pdb_filename)!=EOF) {

	    strcpy(pdb_filename,next_pdb_filename);
	    strcat(pdb_filename,".pdb");
	    fprintf(stderr,"Reading PDB from %s\n",pdb_filename);
	    if(NULL == freopen(pdb_filename, "r", stdin)){
		stop("Error, PDB file cannot be opened\n");
	    }

	    /* CREATE PDB CHAIN (_CHAINT) AND BIASMAP LIBRARY */
	    int i = 1;
	    for (i = 1; pdbin(temporary,&sim_params,stdin) != EOF; i++) { // all should only have 1 chain in the PDBs!!

		if (i==2) stop("All PDBs should have only one chain.");
		n_proteins++;

		/* allocate memory */
		/* library */
		all_chains = (Chain*)realloc(all_chains,n_proteins * sizeof(Chain));
		all_chains[n_proteins-1].aa=NULL; all_chains[n_proteins-1].xaa=NULL; all_chains[n_proteins-1].erg=NULL; all_chains[n_proteins-1].xaa_prev=NULL;
		allocmem_chain(&(all_chains[n_proteins-1]), temporary->NAA, temporary->Nchains);
		/* chaint-s for simulate */
		all_chaints = (Chaint*)realloc(all_chaints,n_proteins * sizeof(Chaint));
		all_chaints[n_proteins-1].aat = NULL; all_chaints[n_proteins-1].xaat = NULL; all_chaints[n_proteins-1].ergt = NULL; all_chaints[n_proteins-1].xaat_prev = NULL; 

		/* copy temporary into the main one */
		copybetween(&((all_chains)[n_proteins-1]),temporary);
		freemem_chain(temporary);

		/* project the peptide onto the CRANKITE model */
		if(sim_params.protein_model.fixit) fixpeptide((all_chains)[n_proteins-1].aa, (all_chains)[n_proteins-1].NAA, &(sim_params.protein_model));
		chkpeptide((all_chains)[n_proteins-1].aa, (all_chains)[n_proteins-1].NAA, &(sim_params.protein_model));
		initialize(&(all_chains[n_proteins-1]),&(all_chaints[n_proteins-1]),&sim_params); // peptide modification

		/* contact map library */
		all_biasmaps = (Biasmap*)realloc(all_biasmaps,n_proteins*sizeof(Biasmap));
		all_biasmaps[n_proteins-1].distb = NULL;
		sim_params.protein_model.contact_map_file = realloc(sim_params.protein_model.contact_map_file,DEFAULT_LONG_STRING_LENGTH*sizeof(char));
		strcpy(sim_params.protein_model.contact_map_file,next_pdb_filename);
		strcat(sim_params.protein_model.contact_map_file,".icm");
		fprintf(stderr,"Reading in contact map from file %s\n",sim_params.protein_model.contact_map_file);
		biasmap_initialise(&(all_chains[n_proteins-1]),&(all_biasmaps[n_proteins-1]),&(sim_params.protein_model));
		//energy matrix init must be called once all parameters have been fixed
		//energy_matrix_calculate(&(all_chains[n_proteins-1]),&(all_biasmaps[n_proteins-1]),&(sim_params.protein_model));

	    }
	    if (i==1) stop("ERROR! EOF while reading in from input PDB file.");

	}
	free(temporary);
	fclose(list_file);
	list_file = NULL;

	/* the biasmaps are different for the Gly residues (no beta-contacts) */
	//for (i=0; i<n_proteins; i++) {
	//	fprintf(stderr,"Biasmap of protein %d:\n",i);
	//	print_contact_map(all_biasmaps[i].distb,all_biasmaps[i].NAA,stderr);
	//}


	/* NEED A LIBRARY OF SIMULATION_PARAMS AND A LIBRARY COPY TO DO THE MC-S */
	fprintf(stderr,"Initialising simulation parameters...\n");
	/* copy sim_params */
	simulation_params *sim_params_sim = NULL;
	sim_params_sim = (simulation_params*)realloc(sim_params_sim,n_proteins * sizeof(simulation_params));
	for (int i=0; i<n_proteins; i++) {
	   sim_params_copy(&((sim_params_sim)[i]),&sim_params);
	}
	param_finalise(&sim_params);

	fprintf(stderr,"Initialising library copy...\n");
	/* initialise, but not yet copy of the library for simulate */
	all_chains_sim = (Chain*)realloc(all_chains_sim,n_proteins * sizeof(Chain));
	for (int i=0; i<n_proteins; i++) {
	   all_chains_sim[i].aa=NULL; all_chains_sim[i].xaa=NULL; all_chains_sim[i].erg=NULL; all_chains_sim[i].xaa_prev=NULL;
	   allocmem_chain(&(all_chains_sim[i]), all_chains[i].NAA, all_chains[i].Nchains);
	}


	/* CD LEARN INITIALISATION */
	fprintf(stderr,"Initialising CD learn parameters...\n");
	cd_learn_param_initialise(&cd_params);
	double *energy_change = NULL;
#ifdef DEBUG
//	double total_energy_change;
#endif
	energy_change = (double *)realloc(energy_change,n_proteins * sizeof(double));

	time_t timer2 = time(NULL);
	fprintf(stderr,"Start-up time: %g\n",(double)timer2-(double)timer1);

	/* INITIALISE OUTPUT FILE */
	FILE *output;
	if (cd_params.output_filename == NULL) {
		cd_params.output_filename = realloc(cd_params.output_filename,DEFAULT_LONG_STRING_LENGTH*sizeof(char));
		strcpy(cd_params.output_filename,"out.");
		strcat(cd_params.output_filename,cd_params.pdb_list_filename);
	}
	fprintf(stderr,"Writing output into %s\n",cd_params.output_filename);
	output = fopen(cd_params.output_filename,"a");
	/* print argument line */
	fprintf(output,"#");
	for (int i=0; i<argc; i++) {
		fprintf(output,"%s ",argv[i]);
	}
	fprintf(output,"\n");
	//fprintf(output,"#crankite/cdlearn -l %s -L %s -o %s -K %d -t %x -i %d -p %s -R %s -a %g\n",
	//	cd_params.learn_string,
	//	cd_params.pdb_list_filename,
	//	cd_params.output_filename,
	//	sim_params_sim[0].pace,
	//	sim_params_sim[0].tmask,
	//	cd_params.iter_max,
	//	sim_params_sim[0].prm,
	//	cd_params.restart_filename,
	//	cd_params.adaptive_learning_param);
	fprintf(stderr,"The initial parameters are: ");
	for (int i=0; i<36; i++) {
	    fprintf(output,"%g ",cd_params.aa[i]);
	    fprintf(stderr,"%g ",cd_params.aa[i]);
	}
	for (int i=0; i<36; i++) {
	    //fprintf(output,"%g ",cd_params.raa[i]);
	    fprintf(output,"%g ",cd_params.total_energy_gradient[i]);
	    //fprintf(stderr,"%g ",cd_params.raa[i]);
	    fprintf(stderr,"%g ",cd_params.total_energy_gradient[i]);
	}
	fprintf(output,"\n");


	/* CD LEARN ITERATIONS */
	fprintf(stderr,"CD learning iterations...\n");
	char next_restart_filename[DEFAULT_LONG_STRING_LENGTH];
	char index[10];
	if (cd_params.restart_filename==NULL) {
		copy_string(&(cd_params.restart_filename),"cdlearn.restart");
	}
	for (int iter = 0; iter < cd_params.iter_max+1; iter++) {

	    /* print restart file */
	    if (iter%100==0) {
		strcpy(next_restart_filename,cd_params.restart_filename);
		sprintf(index,"_%d",iter+cd_params.iter_start);
		strcat(next_restart_filename,index);
		cd_learn_write_restart_file(&cd_params,next_restart_filename);
		fprintf(stderr,"Printing restart file into %s\n",next_restart_filename);
	    }
	    if (iter==cd_params.iter_max) break;
	    fprintf(stderr,"Iteration %d...\n",iter);

	    /* initialise next step by zeroing global counters */
	   /* calc energy test for every single chain */
	   for (int j=0; j<n_proteins; j++) {
	      update_sim_params_from_cd_learn(&(sim_params_sim[j]),&cd_params);
	   }

	   /* loop over all proteins in the library */
	   #pragma omp parallel for schedule(guided,1)
	   for (int j=0; j<n_proteins; j++) {
	      /* 1. update the sim parameters */
	      /*    zero the energy probe parameters before the MC run */
	      for (int k=0; k<36; k++){
		sim_params_sim[j].energy_gradient[k] = 0;
		sim_params_sim[j].energy_probe_1_this[k] = 0;
		sim_params_sim[j].energy_probe_1_last[k] = 0;
	      }
	      /* 2. copy the initial PDB from the original library, incl. energy matrix */
	      copybetween(&(all_chains_sim[j]),&((all_chains)[j]));
	      /*    energy matrix must agree with the parameters */
	      initialize_sidechain_properties(&(sim_params_sim[j].protein_model));
	      energy_matrix_calculate(&(all_chains_sim[j]),&(all_biasmaps[j]),&(sim_params_sim[j].protein_model));
	      /* 3. temporary chain is OK, its energy will be recalculated anyway */
	      /* 4. biasmap is OK. */
	      //print_contact_map(all_biasmaps[j].distb,all_biasmaps[j].NAA,output);

	      energy_change[j] = simulate(&(all_chains_sim[j]),&(all_chaints[j]),&(all_biasmaps[j]),&(sim_params_sim[j]),j,iter);

#ifdef DEBUG
//	      fprintf(stderr,"The energy gradient for protein %d is: ",j);
//	      for (int k=0; k<36; k++) {
//		fprintf(stderr,"%g ",sim_params_sim[j].energy_gradient[k]);
//	      }
//	      fprintf(stderr,"\n");
//	      fprintf(stderr,"The energy gradient calculation for protein %d is: \n",j);
//	      for (int k=0; k<36; k++) {
//		fprintf(stderr,"%g = %g - %g\n",sim_params_sim[j].energy_gradient[k], sim_params_sim[j].energy_probe_1_this[k], sim_params_sim[j].energy_probe_1_last[k]);
//	      }
#endif
	   }

	   /* calc gradient of energy for maximum likelihood */
	   for (int k=0; k<36; k++) {
	      cd_params.total_energy_gradient[k] = 0;
	   }
	   for (int k=0; k<36; k++) {
	      /* accumulate the gradient from every chain */
	      for (int j=0; j<n_proteins; j++) {
		 cd_params.total_energy_gradient[k] += sim_params_sim[j].energy_gradient[k];
	      }
	      cd_params.total_energy_gradient[k] /= n_proteins;
	   }

	    /* calculate new values */
	    for (int k=0; k<36; k++) {
		if (cd_params.raa[k] != 0) {
		    if (cd_params.adaptive_learning_param > 0) { /* adaptive ML step */
			fprintf(stderr,"Adaptive learning: %g of the previous move will be reused.\n",cd_params.adaptive_learning_param);
			double move;
			move = cd_params.adaptive_learning_param * cd_params.previous_move[k] + cd_params.total_energy_gradient[k] * cd_params.raa[k];
			cd_params.aa[k] += move;
			cd_params.previous_move[k] = move;
		    } else { /* standard ML step */
			cd_params.aa[k] += cd_params.total_energy_gradient[k] * cd_params.raa[k];
		    }
		}
	    }
	    /* check that new values are reasonable (e.g. VDW radii are positive), if not, cry and adjust them */
	    correct_negative_aa_values(&cd_params);

	    /* print into output file*/
	    for (int k=0; k<36; k++) { /* parameters */
		fprintf(output,"%g ",cd_params.aa[k]);
	    }
	    for (int k=0; k<36; k++) { /* energy gradient wrt parameters */
		fprintf(output,"%g ",cd_params.total_energy_gradient[k]);
	    }
	    fprintf(output,"\n");

#ifdef DEBUG
//	    /* out of curiosity, print total energy change */
//	    total_energy_change = 0;
//	    for (int j=0; j<n_proteins; j++) {
//		total_energy_change += abs(energy_change[j]);
//	    }
//	    fprintf(stderr,"Total energy change of all proteins %g\n",total_energy_change);
#endif

	    /* flush the output */
	    fflush(output);

	    time_t timer3 = time(NULL);
	    fprintf(stderr," %g s\n",(double)timer3-(double)timer2);
	    timer2 = timer3;

	}

	fclose(output);
	output = NULL;
	fprintf(stderr,"FINISHED!");
	/* finalise memory */
	for (int i=0; i<n_proteins; i++) {

		freemem_chain(&(all_chains[i]));
		freemem_chain(&(all_chains_sim[i]));
		freemem_chaint(&(all_chaints[i]));

		free(all_biasmaps[i].distb);

		param_finalise(&sim_params_sim[i]);
	}
	free(all_chains);
	free(all_chains_sim);
	free(all_chaints);
	free(all_biasmaps);
	free(sim_params_sim);
	free(energy_change);
	cd_param_finalise(&cd_params);

	return EXIT_SUCCESS;
}
