/*
 * flex.h
 *
 *  Created on: 14 Sep 2012
 *      Author: nik
 */

#ifndef FLEX_H_
#define FLEX_H_

#ifdef PARALLEL
#include<mpi.h>
void ns_for_flex_processor(MPI_Comm FLEX_WORLD,int rank,Biasmap *biasmap,simulation_params* sim_params);
void update_flex_amplitude(int flex_iteration,simulation_params * flex_params,FLEX_data * flex_data);
int check_finished(simulation_params *sim_params);
#endif

void initialize_flex(Chain *chain,Chain**input_chains,Biasmap *biasmap,simulation_params* sim_params,double **rmsd);
int output_and_run_flex(Chain *chain,Biasmap *biasmap,Chain *input_chains,simulation_params* sim_params,double *rmsd);
int read_in_after_flex(Chain *chain,Biasmap *biasmap,Chain *input_chains,simulation_params *sim_params,double *rmsd);
void finalize_flex(Chain *chain, Chain **input_chains,simulation_params *sim_params,double **rmsd);


int check_flex_ready(simulation_params *sim_params);


void output_covin(Chain *chain, model_params *protein_model,char *filename);
void output_stackedin(char* filename);
void Hyd_pdbout(Chain *chain,model_params *protein_model,char *filename,char *ca_filename);
void output_hphobesin(Chain * chain, model_params *protein_model,char *filename, FILE *logfile);
void output_hbondsin(Chain * chain, model_params *protein_model,int use_bias_only,char *filename,double nma_hstrength, FILE *logfile);


#endif /* NMA_H_ */
