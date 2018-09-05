/*
 * nma.h
 *
 *  Created on: 14 Sep 2012
 *      Author: nik
 */

#ifndef NMA_H_
#define NMA_H_

#ifdef PARALLEL
#include<mpi.h>
void ns_for_nma_processor(MPI_Comm NMA_WORLD,int rank,Biasmap *biasmap,simulation_params* sim_params);
#endif

void initialize_nma(Chain *chain,Chain**input_chains,simulation_params* sim_params);
int output_and_run_nma(Chain *chain,Biasmap *biasmap,Chain *input_chains,simulation_params* sim_params);
int read_in_after_nma(Chain *chain,Biasmap *biasmap,Chain *input_chains,simulation_params *sim_params);
void finalize_nma(Chain *chain, Chain **input_chains,simulation_params *sim_params);


int check_NMA_ready(simulation_params *sim_params);


void output_covin(Chain *chain, model_params *protein_model,char *filename);
void output_stackedin(char* filename);
void Hyd_pdbout(Chain *chain,model_params *protein_model,char *filename,char *ca_filename);
void output_hphobesin(Chain * chain, model_params *protein_model,char *filename);
void output_hbondsin(Chain * chain, model_params *protein_model,char *filename,double nma_hstrength);


#endif /* NMA_H_ */
