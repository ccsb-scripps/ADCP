/* Nested Sampling procedure for Crankite
 * 
 * Skeleton of code from J.Skilling p188 of
 * Data Analysis A Bayesian Tutorial D.S.Sivia
 *
 * Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
 * 
 */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<time.h>
#include"error.h"
#include"params.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"
#include"vdw.h"
#include"energy.h"
#include"metropolis.h"
#include"probe.h" 
#include"nested.h"

#ifdef PARALLEL
#include<mpi.h>
#include"random16.h"
#include"flex.h"
#endif
#include"checkpoint_io.h"

int lowtemp = 0;

#define PLUS(x,y) (x>y ? x + log(1+exp(y-x)) : y + log(1+exp(x-y)))


/* A bit of ordering of the chains wrt. the logL in ascending order.
   It takes the first one, and moves it towards the end a bit, or all
   the way to the end, if it has the largest logL.
   If the first one has the largest logL, it will be moved to the end. */
void heapifyhashout(ChainHash* chainhash ,int length){
  int booly = 0;
  int k = 1,i;
  ChainHash v;
  copyhash(&v,&(chainhash[k]));
  while((booly == 0) && 2 * k <= length){
    i = 2 * k;
    if(i < length) //2 children
      if(chainhash[i].ll > chainhash[i+1].ll) i += 1;
    if(v.ll < chainhash[i].ll) booly = 1;
    else{copyhash(&(chainhash[k]),&(chainhash[i])); k = i;} 
  }	  
  copyhash(&(chainhash[k]),&v);	
} 


/* A bit of ordering of the chains wrt. the logL in ascending order.
   It takes the last one that has probably just been added to the heap,
   and moves it down the line a bit, or to the front, if it has the
   smallest logL.
   If the last one has the smallest logL, it will be moved to the front. */
void heapifyhashin(ChainHash *chainhash, int length){
  int current = length;
  int booly = 0;
  ChainHash v;
  while((booly == 0) && current > 1){
    int next = current / 2;
    if(chainhash[current].ll < chainhash[next].ll){
      copyhash(&v,&(chainhash[next]));
      copyhash(&(chainhash[next]),&(chainhash[current]));
      copyhash(&(chainhash[current]),&v);
      current /= 2;
    }
    else{
      booly = 1;
    }
  }	
}

void repairheap(ChainHash* chainhash,int N){
  constructhashheap(chainhash, N);
}


void new_amplitude(Chain **cpoints, Biasmap *biasmap, int current_stored, simulation_params *sim_params, void *mpi_comm){

  int rank = 0;
  int P = 1;
  double amptot = 0;
  double oldamp = sim_params->amplitude;
  Chaint *chaint;
  chaint = (Chaint*)malloc(sizeof(Chaint));
  chaint->aat = NULL; chaint->ergt = NULL;  chaint->xaat = NULL; chaint->xaat_prev = NULL;
  double currE;
  int copies;
  Chain temporary;
  temporary.aa = NULL; temporary.xaa = NULL; temporary.erg = NULL; temporary.xaa_prev = NULL;
  allocmem_chain(&temporary,(*cpoints)[0].NAA,(*cpoints)[0].Nchains);

  //calculate amplitude
#ifdef PARALLEL
  MPI_Comm *MPI_COMM = mpi_comm;
  MPI_Status status;
  MPI_Comm_size(*MPI_COMM, &P);
  MPI_Comm_rank(*MPI_COMM, &rank);
  if(P == 1){
    for(int k = 0; k < 5; k++){
      copies = (int)(rand()/(double)RAND_MAX * current_stored) % (current_stored);
      sim_params->amplitude = oldamp;
      copybetween(&temporary,&((*cpoints)[copies]));
      aat_init(&temporary,chaint);
      currE = -temporary.ll;
      move(&temporary,chaint,biasmap,sim_params->logLstar,&currE,-1, sim_params);
      for(int i = 1; (i < sim_params->iter_max || i < 1024); i++){
        move(&temporary,chaint,biasmap,sim_params->logLstar,&currE,1, sim_params);
      }
      amptot += sim_params->amplitude;
    }
    sim_params->amplitude = amptot / 5;
    return;
  }


  if (rank != 0) { // one chain per processor
#else
    for(int k = 0; k < 5; k++){ //5 chains only
#endif
      copies = (int)(rand()/(double)RAND_MAX * current_stored) % (current_stored);

      //recalculate amplitude on temporary, not changing the actual NS points
      copybetween(&temporary,&((*cpoints)[copies]));
      aat_init(&temporary,chaint);
#ifndef PARALLEL
      sim_params->amplitude = oldamp;
#endif
      currE = -temporary.ll;
      move(&temporary,chaint,biasmap,sim_params->logLstar,&currE,-1, sim_params);
      for(int i = 1; (i < sim_params->iter_max || i < 1024); i++){
        move(&temporary,chaint,biasmap,sim_params->logLstar,&currE,1, sim_params);
      }

      //collect and add up amplitude
#ifdef PARALLEL
      //send my new amplitude
      MPI_Send(&(sim_params->amplitude),1,MPI_DOUBLE,0,sim_params->iter*6+5,*MPI_COMM);
    } else {
      for(int k = 1; k < P;k++){
        MPI_Recv(&(sim_params->amplitude),1,MPI_DOUBLE,k,sim_params->iter*6+5,*MPI_COMM,&status);
        amptot += sim_params->amplitude;
      }
      sim_params->amplitude = amptot / (P-1);
    }
#else
    amptot += sim_params->amplitude;
  }
  sim_params->amplitude = amptot / 5;
  //fprintf(stderr,"Amplitude: %f\n",-sim_params->amplitude);
#endif

  //broadcast amplitude
#ifdef PARALLEL
  //receive the new amplitude from rank 0
  MPI_Bcast(&(sim_params->amplitude),1,MPI_DOUBLE,0,*MPI_COMM);
#endif

  freemem_chain(&temporary);
  freemem_chaint(chaint);
  free(chaint);
}

#ifdef PARALLEL
/*set up 2 separate communicators, FLEX_WORLD for
 * processors 0,1,...nma_params.number_of_processors >> the master processor and the
 * (for now) 1 processor involved in FLEX
 * and NS_WORLD for processors
 * 0,nma_params.number_of_processors+1,FLEX_params.number_of_processors+2,...P-1
 * for the standard NS part
 */
void setup_communicators(MPI_Comm * FLEX_WORLD, MPI_Comm * NS_WORLD, int S, int *P, int *rank,int *SP, int *Srank, int *NSgrouprank, int *FLEXgrouprank, int *in_NS){
  MPI_Comm_size(MPI_COMM_WORLD, P);
  MPI_Comm_rank(MPI_COMM_WORLD, rank);

  if(*P-S < 1){
    stop("You must have at least 1 processor running the NS simulation\n");
  }

  int setup_groups;
  MPI_Group NS_GROUP, FLEX_GROUP, ALL_GROUP;

  MPI_Comm_group( MPI_COMM_WORLD, &ALL_GROUP );

  if(S > 0){
    int *FLEX_group_ranks = malloc(sizeof(int)*(S+1));

    for(setup_groups = 0; setup_groups <= S; setup_groups++ ){
      FLEX_group_ranks[setup_groups] = setup_groups;
    }
    MPI_Group_incl(ALL_GROUP,S+1,FLEX_group_ranks,&FLEX_GROUP);
    MPI_Comm_create(MPI_COMM_WORLD,FLEX_GROUP,FLEX_WORLD);
    free(FLEX_group_ranks);
  }

  int *NS_group_ranks = malloc(sizeof(int)*(*P-S));
  NS_group_ranks[0] = 0;
  for(setup_groups = 1; setup_groups < *P-S; setup_groups++ ){
    NS_group_ranks[setup_groups] = setup_groups+S;
  }
  MPI_Group_incl(ALL_GROUP,*P-S,NS_group_ranks,&NS_GROUP);
  MPI_Comm_create(MPI_COMM_WORLD,NS_GROUP,NS_WORLD);
  free(NS_group_ranks);

  if(S>0 && *rank <= S){
    MPI_Group_rank(FLEX_GROUP,FLEXgrouprank);
    MPI_Comm_size(*FLEX_WORLD, SP);
    MPI_Comm_rank(*FLEX_WORLD, Srank);
    if(*Srank != 0) *in_NS = -1;
  }
  if(*rank == 0 || *rank > S){
    MPI_Group_rank(NS_GROUP,NSgrouprank);
    MPI_Comm_size(*NS_WORLD, P);
    MPI_Comm_rank(*NS_WORLD, rank);
  }

}

void deal_with_flex(int N, int finished,int rank, simulation_params *sim_params, MPI_Comm NS_WORLD, MPI_Comm FLEX_WORLD, ChainHash *chainhash, Chain* temporary,Chain *cpoints){
  MPI_Status status;


  int flex_status_for_ns;
  int one = 1;
  if(rank == 0){


    int flex_ready = check_flex_ready(sim_params);
    if(flex_ready == 1 || sim_params->iter == sim_params->iter_start){
      if(flex_ready == 1){

        flex_status_for_ns = 1;
        MPI_Bcast(&flex_status_for_ns,1,MPI_INT,0,NS_WORLD);
        MPI_Send(&(sim_params->logLstar),1,MPI_DOUBLE,one,10,FLEX_WORLD);


        int number_to_recv = 0;
        MPI_Recv(&number_to_recv,1,MPI_INT,one,11,FLEX_WORLD,&status);
        MPI_Bcast(&number_to_recv,1,MPI_INT,0,NS_WORLD);

        int receives;

        for(receives = 0; receives < number_to_recv; receives++){

          mpi_rec_chain(temporary, one, 0, &sim_params->logLstar, 0, FLEX_WORLD);
          int overwrite;
          overwrite = 1 + ((int)(rand()/(double)RAND_MAX * (N))) % (N);


          chainhash[overwrite].ll = temporary->ll;


          int where_to_send = chainhash[overwrite].processor;


          MPI_Bcast(&where_to_send,1,MPI_INT,0,NS_WORLD);
          if(where_to_send == 0){
            copybetween(&cpoints[chainhash[overwrite].index],temporary);

          }
          else{
            MPI_Send(&(chainhash[overwrite].index),1,MPI_INT, chainhash[overwrite].processor ,sim_params->iter,NS_WORLD);
            mpi_send_chain(temporary,0,where_to_send,&sim_params->logLstar,sim_params->iter,NS_WORLD);
          }

        }
        repairheap(chainhash,N);



      }

      else{
        flex_status_for_ns = 0;
        MPI_Bcast(&flex_status_for_ns,1,MPI_INT,0,NS_WORLD);
      }

      int tosend = 1;//1 + ((int)(rand()/(double)RAND_MAX * (N))) % (N);

      int where_from = chainhash[tosend].processor;


      MPI_Bcast(&where_from,1,MPI_INT,0,NS_WORLD);
      if(where_from == 0){
        copybetween(temporary,&cpoints[chainhash[tosend].index]);
      }
      else{
        MPI_Send(&(chainhash[tosend].index),1,MPI_INT, chainhash[tosend].processor ,sim_params->iter,NS_WORLD);
        mpi_rec_chain(temporary, chainhash[tosend].processor, 0, &sim_params->logLstar, sim_params->iter, NS_WORLD);
      }
      //chain sent to FLEX
      mpi_send_chain(temporary,0,one,&sim_params->logLstar,1,FLEX_WORLD);
    }
    else{
      flex_status_for_ns = -1; //done
      MPI_Bcast(&flex_status_for_ns,1,MPI_INT,0,NS_WORLD);
    }



  }
  else{ //not rank 0
    MPI_Bcast(&flex_status_for_ns,1,MPI_INT,0,NS_WORLD);

    if(flex_status_for_ns == 1){
      int number_to_recv = 0;
      MPI_Bcast(&number_to_recv,1,MPI_INT,0,NS_WORLD);

      int receives;
      for(receives = 0; receives < number_to_recv; receives++){
        int want_it;
        MPI_Bcast(&want_it,1,MPI_INT,0,NS_WORLD);
        if(want_it == rank){
          int which_index;
          MPI_Recv(&(which_index),1,MPI_INT,0,sim_params->iter,NS_WORLD,&status);
          mpi_rec_chain(&(cpoints[which_index]), 0, rank, &sim_params->logLstar , sim_params->iter,NS_WORLD);
        }
      }
    }
    if(flex_status_for_ns >= 0){
      int where_from;
      MPI_Bcast(&where_from,1,MPI_INT,0,NS_WORLD);


      if(where_from == rank){
        int which_index;
        MPI_Recv(&which_index,1,MPI_INT, 0 ,sim_params->iter,NS_WORLD,&status);
        mpi_send_chain(&cpoints[which_index],rank,0,&sim_params->logLstar,sim_params->iter,NS_WORLD);
      }
    }



  }
  if(finished == 0 && rank == 0){

    char filename[DEFAULT_SHORT_STRING_LENGTH];
    //strcpy(filename,sim_params->flex_params.output_path);
    strcpy(filename,"finished");
    FILE *fptr = fopen(filename,"w");
    fclose(fptr);

  }

}


#endif


void check_to_output_checkpoint_file(Chain* cpoints, int current_stored, int N, simulation_params *sim_params, int rank, void*comm_pointer){
  //if outputting checkpoint file
  // The output printing has to happen before updating the logX, weight, logZ etc., otherwise
  //   a) the logZ will contain the highest energy configuration's weight twice
  //   b) logX will have been advanced using the old run's alpha and we cannot correct for it any more (matters if P changes)
  // This means we print the previous logX, logZ, logL*, H and amplitude
  if(sim_params->checkpoint == 1){
    sim_params->N = N;
    if(sim_params->iter != sim_params->iter_start) {
      //fprintf(stderr,"trying to print checkpoint\n");
      output_checkpoint_file(cpoints, current_stored, sim_params,comm_pointer);
    }
    if (rank == 0) {
      if( sim_params->iter % sim_params->num_NS_per_checkpoint == 0 ){
        char *out = (char*)malloc(sizeof(char)*1010);
        if(sim_params->outfile_name != NULL){
          sprintf(out,"%s_%d",sim_params->outfile_name,sim_params->checkpoint_counter);
          if (sim_params->outfile) fclose(sim_params->outfile);
          sim_params->outfile = fopen(out, "w");
        }
        free(out);
      }
    }
  }
}

void find_worst(simulation_params *sim_params, ChainHash *chainhash, int *heaplength, int N, int P ){
  /*find P worst samples in positions, N,N-1,N-2,... */
  for(int k = 0; k < P; k++){
    copyhash(&(chainhash[0]),&(chainhash[1]));
    copyhash(&(chainhash[1]),&(chainhash[*heaplength]));
    copyhash(&(chainhash[*heaplength]),&(chainhash[0]));
    (*heaplength)--;
    heapifyhashout(chainhash,*heaplength);
  }

  sim_params->logLstar = chainhash[N-P+1].ll;
}


void update_NS_parameters(simulation_params *sim_params, int *converged, double *logZnew, double *lweight){
  /*update weight = width * likelihood */
  *lweight = sim_params->log_DeltaX + sim_params->logLstar*sim_params->thermobeta;
  //fprintf(stderr,"LWJ %f %d\n",sim_params->log_DeltaX,sim_params->iter);
  sim_params->logX = sim_params->logX_start + sim_params->Delta_logX*(sim_params->iter-sim_params->iter_start+1);

  /* update evidence and information */
  *logZnew = PLUS(sim_params->logZ,*lweight);
  sim_params->H = exp(*lweight - *logZnew) * sim_params->logLstar * sim_params->thermobeta  + exp(sim_params->logZ-*logZnew)*(sim_params->H+sim_params->logZ) - *logZnew;

  /* control */
  if(sim_params->checkpoint == 1 && sim_params->iter % sim_params->num_NS_per_checkpoint == 0 && sim_params->iter != sim_params->iter_start) {
    fprintf(stderr,"%d %lf %lf %lf %lf %lf\n", sim_params->iter,    sim_params->logX,   sim_params->logLstar,sim_params->logZ,sim_params->H,sim_params->amplitude);
  }

  /*Check if converged */
  if(sim_params->lowtemp == 1 && fabs(*logZnew - sim_params->logZ) < 1e-8)  *converged = 0;
  sim_params->logZ = *logZnew;
}

#ifdef PARALLEL

void MC_first(ChainHash *ChainHash,Chain *cpoints, Chaint* chaint, int current_stored, Biasmap *biasmap,simulation_params *sim_params, int rank, int N,int P, MPI_Comm *NSWORLD){
  int i;
  MPI_Status status;


  double *newls = (double*)malloc(sizeof(double)*(1+N/P));
  for(i = 0; i < current_stored; i++){
    int iter;
    double currE = -cpoints[i].ll;
    for(iter= 0; iter < sim_params->number_initial_MC; iter++){
      move(&(cpoints[i]),chaint,biasmap,0,&currE,0,sim_params);
    }
    cpoints[i].ll = -currE;
    newls[i] = cpoints[i].ll;
  }

  if(rank == 0){
    fprintf(stderr,"Beware: %d initial MC moves before NS started\n",sim_params->number_initial_MC);
    for(i = 0; i < current_stored; i++){
      ChainHash[1+i*P].ll = newls[i];
    }
    int proc;
    if(P > 1){
      for(proc = 1; proc < P; proc++){
        MPI_Recv(newls,1+N/P,MPI_DOUBLE,proc,proc,*NSWORLD,&status);
        for(i = 0; i < N/P+1; i++){
          if(1+i*P+proc <= N) ChainHash[1+i*P+proc].ll = newls[i];
        }
      }

    }
  }
  else{
    MPI_Send(newls,1+N/P,MPI_DOUBLE,0,rank,*NSWORLD);
  }

  free(newls);

}

#ifndef FAST
void collect_chains(ChainHash *chainhash,Chain *cpoints, Chain*chaincopies, simulation_params *sim_params,  int rank, int P, int N,MPI_Comm *NS_WORLD,Instructions *instructions){
  //first collect the P chains which will be the start of the MMC


  int copies, got_all;
  MPI_Status status;
  int minus1 = -1;

  //collect the chains

  if (rank == 0) {
    //random number generation would fail if N==P
    if (N==P) stop("number of processors = number of active points.  All chains would be scrapped at each MC step.");
    if (N<P) stop("number of processors > number of active points.  All chains would be scrapped at each MC step.");
    //first collect the P chains which will be the start of the MMC
    int * start_chains = (int*)malloc(sizeof(int)*P);
    for(int k = 0; k < P; k++){
      copies = 1 + ((int)(rand()/(double)RAND_MAX * (N-P))) % (N-P);
      start_chains[k] = copies;
    }
    for(int k = 0; k < P; k++){

      if(chainhash[start_chains[k]].processor == 0){
        copybetween(&(chaincopies[k]),&(cpoints[chainhash[start_chains[k]].index]));
      }
      else{
        MPI_Send(&(chainhash[start_chains[k]].index),1,MPI_INT,chainhash[start_chains[k]].processor,sim_params->iter,*NS_WORLD);
        mpi_rec_chain(&(chaincopies[k]), chainhash[start_chains[k]].processor, 0,&sim_params->logLstar , sim_params->iter,*NS_WORLD);
      }


    }
    //let all processors know we're done
    for(int i = 1; i < P; i++){
      MPI_Send(&minus1,1,MPI_INT,i,sim_params->iter,*NS_WORLD);
    }

    //then send one to each processor
    for(int k = 1; k < P; k++){
      mpi_send_chain(&(chaincopies[k]), 0, k,&sim_params->logLstar,sim_params->iter,*NS_WORLD);
    }
    free(start_chains);
  } //end rank 0

  else{

    //processor 0 is asking for a random  P cpoints spread throughout all processors
    //send any it wants
    got_all = 0;
    while(got_all != -1){
      MPI_Recv(&got_all,1,MPI_INT,0,sim_params->iter,*NS_WORLD,&status);
      if(got_all != -1)   mpi_send_chain(&(cpoints[got_all]), rank,0, &sim_params->logLstar , sim_params->iter,*NS_WORLD);
    }
    //receive the chain I'm expected to MC
    mpi_rec_chain(&chaincopies[0],0,rank,&sim_params->logLstar,sim_params->iter,*NS_WORLD);

  }
}

void return_and_reheap_chains(ChainHash *chainhash,Chain *cpoints, Chain*chaincopies, simulation_params *sim_params,  int rank, int P, int *heaplength, MPI_Comm *NS_WORLD, Instructions *instructions){
  //recieve the P chains on the master processor
  int minus1 = -1;
  int got_all;
  MPI_Status status;
  if (rank == 0) {
    for(int k = 1; k < P; k++){
      mpi_rec_chain(&(chaincopies[k]), k, 0,&sim_params->logLstar,sim_params->iter,*NS_WORLD);
    }

  }
  else {
    mpi_send_chain(&chaincopies[0],rank,0,&sim_params->logLstar,sim_params->iter,*NS_WORLD);
  }

  if (rank == 0) {
    //put them back in heap and send them to overwrite the worst ones
    for(int k = 0; k < P; k++){
      (*heaplength)++;
      chainhash[*heaplength].ll = chaincopies[k].ll;
      if(chainhash[*heaplength].processor == 0){
        copybetween(&cpoints[chainhash[*heaplength].index],&chaincopies[k]);
      }
      else{
        MPI_Send(&(chainhash[*heaplength].index),1,MPI_INT, chainhash[*heaplength].processor ,sim_params->iter,*NS_WORLD);
        mpi_send_chain(&(chaincopies[k]),0, chainhash[*heaplength].processor,&sim_params->logLstar,sim_params->iter,*NS_WORLD);
      }
      heapifyhashin(chainhash,*heaplength);
    }

    //tell processors we've given out all chains
    for(int i = 1; i < P; i++){
      MPI_Send(&minus1,1,MPI_INT,i,sim_params->iter,*NS_WORLD);
    }
  }
  else {


    //get given any of the new points to replace the worst P points in the collective cpoints arrays
    got_all = 0;
    while(got_all != -1){
      MPI_Recv(&got_all,1,MPI_INT, 0 ,sim_params->iter,*NS_WORLD,&status);
      if(got_all != -1) mpi_rec_chain(&(cpoints[got_all]),0,rank,&sim_params->logLstar,sim_params->iter,*NS_WORLD);
    }
  }

}


#else

void initialize_instruction_set(Instructions * instruction, int P){
  instruction->length = 6*P+1;
  instruction->current_position = -1;
  instruction->instructions = (int*)malloc(sizeof(int)*instruction->length);
}

void finalize_instruction_set(Instructions * instruction){
  free(instruction->instructions);
}

void carry_out_first_instructions(ChainHash *chainhash,Chain *cpoints, Chain*chaincopies, simulation_params *sim_params,  int rank, int P, int N,MPI_Comm *NS_WORLD,Instructions *instructions){
  instructions[0].current_position = -3;
  int instruct; int proc;
  int index;
  do{
    instructions[0].current_position+=3;
    instruct = instructions[0].instructions[instructions[0].current_position];
    proc = instructions[0].instructions[instructions[0].current_position+1];
    index = instructions[0].instructions[instructions[0].current_position+2];
    if(instruct == 1 ){ //send
      mpi_send_chain(&(cpoints[index]), rank, proc,&sim_params->logLstar,sim_params->iter,*NS_WORLD);
    }
    else if(instruct == -1){ //receive
      mpi_rec_chain(&(chaincopies[0]), proc,rank,&sim_params->logLstar , sim_params->iter,*NS_WORLD);
    }
  }while(instruct != 0);
  if(index != -1){
    copybetween(&(chaincopies[0]),&(cpoints[index]));
  }
  instructions[0].current_position+=3;
}

void carry_out_last_instructions(ChainHash *chainhash,Chain *cpoints, Chain*chaincopies, simulation_params *sim_params,  int rank, int P, int N,MPI_Comm *NS_WORLD,Instructions *instructions){
  instructions[0].current_position-=3;
  int instruct; int proc;
  int index;
  do{
    instructions[0].current_position+=3;
    instruct = instructions[0].instructions[instructions[0].current_position];
    proc = instructions[0].instructions[instructions[0].current_position+1];
    index = instructions[0].instructions[instructions[0].current_position+2];
    if(instruct == 1 ){ //send
      mpi_send_chain(&(chaincopies[0]), rank, proc,&sim_params->logLstar,sim_params->iter,*NS_WORLD);
    }
    else if(instruct == -1){ //receive
      mpi_rec_chain(&(cpoints[index]), proc,rank,&sim_params->logLstar , sim_params->iter,*NS_WORLD);
    }
  }while(instruct != 0);
  if(index != -1){
    copybetween(&(cpoints[index]),&(chaincopies[0]));
  }
}

void output_instructions(Instructions *instructions, int rank){
  fprintf(stderr, "Rank %d Number Instructions %d:\n",rank,instructions->current_position/3);
  for(int i = 0; i < instructions->current_position/3; i++){
    fprintf(stderr, "%d %d %d\n",instructions->instructions[3*i],instructions->instructions[3*i+1],instructions->instructions[3*i+2]);
  }
}


void send_and_receive_instructions( Instructions * instructions, int rank, int P, MPI_Comm *NS_WORLD){

  if(rank == 0){
    //output_instructions(&(instructions[0]),0);
    for(int k = 1; k < P; k++){
      //output_instructions(&(instructions[k]),k);
      MPI_Send(&(instructions[k].current_position),1,MPI_INT,k,25,*NS_WORLD);
      MPI_Send(instructions[k].instructions,instructions[k].current_position,MPI_INT,k,50,*NS_WORLD);
    }
  }
  else{
    MPI_Status status;
    MPI_Recv(&(instructions[0].current_position),1,MPI_INT,0,25,*NS_WORLD,&status);
    MPI_Recv(instructions[0].instructions,instructions[0].current_position,MPI_INT,0,50,*NS_WORLD,&status);

  }


}




void collect_chains(ChainHash *chainhash,Chain *cpoints, Chain*chaincopies, simulation_params *sim_params,  int rank, int P, int N,MPI_Comm *NS_WORLD,Instructions *instructions){
  MPI_Bcast(&(sim_params->logLstar),1,MPI_DOUBLE,0,*NS_WORLD);
  if (rank == 0) {
    int *procs = malloc(sizeof(int)*P);
    for(int k = 0; k < P; k++){
      procs[k] = -1;
      instructions[k].current_position = 0;
    }

    int *leftovers = malloc(sizeof(int)*(P*2));
    int number_leftovers = 0;

    int copies;
    for(int k = 0; k < P; k++){
      copies = 1 + ((int)(rand()/(double)RAND_MAX * (N-P))) % (N-P);

      int which_proc = chainhash[copies].processor;
      int which_index = chainhash[copies].index;

      if(procs[which_proc] == -1){
        procs[which_proc] = which_index;
      }
      else{
        leftovers[number_leftovers*2] = which_proc;
        leftovers[number_leftovers*2+1] = which_index;
        number_leftovers++;
      }
    }

//    static int number_chain = 0;
//    number_chain += number_leftovers;

    for(int k = 0; k < P; k++){
      if(procs[k] == -1){
        instructions[k].instructions[instructions[k].current_position] = -1;
        instructions[k].instructions[instructions[k].current_position+1] = leftovers[(number_leftovers-1)*2];
        instructions[k].instructions[instructions[k].current_position+2] = leftovers[(number_leftovers-1)*2+1];
        instructions[k].current_position+=3;
        instructions[leftovers[(number_leftovers-1)*2]].instructions[instructions[leftovers[(number_leftovers-1)*2]].current_position] = +1;
        instructions[leftovers[(number_leftovers-1)*2]].instructions[instructions[leftovers[(number_leftovers-1)*2]].current_position+1] = k;
        instructions[leftovers[(number_leftovers-1)*2]].instructions[instructions[leftovers[(number_leftovers-1)*2]].current_position+2] = leftovers[(number_leftovers-1)*2+1];
        instructions[leftovers[(number_leftovers-1)*2]].current_position+=3;
        number_leftovers--;
      }
    }


    for(int k = 0; k < P; k++){
      instructions[k].instructions[instructions[k].current_position] = 0;
      instructions[k].instructions[instructions[k].current_position+1] = 0;
      instructions[k].instructions[instructions[k].current_position+2] = procs[k];
      instructions[k].current_position+=3;

    }

    //now sort out after MC move instructions
    for(int k = 0; k < P; k++)procs[k] = -1;
    number_leftovers = 0;




    for(int k = 0; k < P; k++){
      copies = N-k;

      int which_proc = chainhash[copies].processor;
      int which_index = chainhash[copies].index;

      if(procs[which_proc] == -1){
        procs[which_proc] = which_index;
      }
      else{
        leftovers[number_leftovers*2] = which_proc;
        leftovers[number_leftovers*2+1] = which_index;
        number_leftovers++;
      }
    }

//    number_chain += number_leftovers;
    //fprintf(stderr,"NUMBER CHAINS %d\n",number_chain);

    for(int k = 0; k < P; k++){
      if(procs[k] == -1){
        instructions[k].instructions[instructions[k].current_position] = 1;
        instructions[k].instructions[instructions[k].current_position+1] = leftovers[(number_leftovers-1)*2];
        instructions[k].instructions[instructions[k].current_position+2] = leftovers[(number_leftovers-1)*2+1];
        instructions[k].current_position+=3;
        instructions[leftovers[(number_leftovers-1)*2]].instructions[instructions[leftovers[(number_leftovers-1)*2]].current_position] = -1;
        instructions[leftovers[(number_leftovers-1)*2]].instructions[instructions[leftovers[(number_leftovers-1)*2]].current_position+1] = k;
        instructions[leftovers[(number_leftovers-1)*2]].instructions[instructions[leftovers[(number_leftovers-1)*2]].current_position+2] = leftovers[(number_leftovers-1)*2+1];
        instructions[leftovers[(number_leftovers-1)*2]].current_position+=3;
        number_leftovers--;
      }
    }


    for(int k = 0; k < P; k++){
      instructions[k].instructions[instructions[k].current_position] = 0;
      instructions[k].instructions[instructions[k].current_position+1] = 0;
      instructions[k].instructions[instructions[k].current_position+2] = procs[k];
      instructions[k].current_position+=3;

    }

    free(procs);
    free(leftovers);



  }


  send_and_receive_instructions(instructions,rank,P,NS_WORLD);
  carry_out_first_instructions(chainhash,cpoints, chaincopies, sim_params, rank,  P,  N,NS_WORLD,instructions);
}

void return_and_reheap_chains(ChainHash *chainhash,Chain *cpoints, Chain*chaincopies, simulation_params *sim_params,  int rank, int P, int *heaplength, MPI_Comm *NS_WORLD, Instructions *instructions){
  carry_out_last_instructions(chainhash,cpoints, chaincopies, sim_params, rank,  P, *heaplength,NS_WORLD,instructions);
  if (rank == 0) {
    MPI_Status status;
    for(int k = 0; k < P; k++){
      (*heaplength)++;
      if(k == 0){
        chainhash[*heaplength].ll = chaincopies[0].ll;
        instructions[k].current_position+=2;
      }
      else{
        instructions[k].current_position--;
        MPI_Recv(&(chainhash[*heaplength].ll),1,MPI_DOUBLE,k,75,*NS_WORLD,&status);
      }

      //fprintf(stderr,"KK %d %d %d\n",k,instructions[k].current_position,instructions[k].instructions[instructions[k].current_position]); fflush(stderr);

      if(instructions[k].instructions[instructions[k].current_position] == -1){
        instructions[k].current_position -=3;
        chainhash[*heaplength].processor = instructions[k].instructions[instructions[k].current_position-1];
      }
      else{
        chainhash[*heaplength].processor = k;
      }

      chainhash[*heaplength].index = instructions[k].instructions[instructions[k].current_position];

      heapifyhashin(chainhash,*heaplength);
    }


  }

  else{
    MPI_Send(&(chaincopies[0].ll),1,MPI_DOUBLE,0,75,*NS_WORLD);
  }
}



#endif


#endif


/* Output NS point (in a serial run), or
   collect the NS point through MPI and output it on the master node. */
void output_NS_point(ChainHash *chainhash, Chain *cpoints, Biasmap *biasmap, simulation_params *sim_params, int N,void*mpi_comm) {
	Chain *chain_to_output = NULL;

	int P = 1;
//get the chain to be output
#ifdef PARALLEL
	MPI_Comm *MPI_COMM = mpi_comm;
	int rank = 0;
	MPI_Status status;
	MPI_Comm_size(*MPI_COMM, &P);
	MPI_Comm_rank(*MPI_COMM, &rank);
	Chain temporary; //if the chain needs collecting from another processor

	if (rank == 0) {
	  temporary.aa = NULL; temporary.xaa = NULL; temporary.erg = NULL; temporary.xaa_prev = NULL;
	  allocmem_chain(&temporary,cpoints[0].NAA,cpoints[0].Nchains);
	  //We need to update AA.id and AA.num, because those are not sent through MPI.
	  for (int i=0; i< sim_params->NAA; i++) {
	    temporary.aa[i].id = sim_params->seq[i];
	    temporary.aa[i].num = i;
	  }
	  //Need to actually get the chain first!

	  MPI_Bcast(&chainhash[N-P+1].processor,1,MPI_INT,0,*MPI_COMM);
	  if(chainhash[N-P+1].processor != 0){
	    MPI_Send(&chainhash[N-P+1].index,1,MPI_INT,chainhash[N-P+1].processor,sim_params->iter,*MPI_COMM );
	    mpi_rec_chain(&temporary, chainhash[N-P+1].processor, 0, &sim_params->logLstar,sim_params->iter,*MPI_COMM);
	    chain_to_output = &temporary;
	  }
	  else{
	    chain_to_output = &(cpoints[chainhash[N-P+1].index]);
	  }
	}
	else {
	  int whichP, whichindex;
	  //if it's my chain processor 0 wants then send it
	  MPI_Bcast(&whichP,1,MPI_INT,0,*MPI_COMM);
	  if(whichP == rank){
	    MPI_Recv(&whichindex,1,MPI_INT,0,sim_params->iter,*MPI_COMM,&status );
	    mpi_send_chain( &(cpoints[whichindex]), rank,0, &sim_params->logLstar,sim_params->iter,*MPI_COMM);
	  }
	}
#else
		chain_to_output = &(cpoints[chainhash[N-P+1].index]);
//	chain_to_output = &(cpoints[chainhash[1].index]);
#endif

//do the tests on it
#ifdef PARALLEL
    if (rank == 0) {
#endif
	tests(chain_to_output,biasmap,sim_params->tmask, sim_params, 0x11, NULL);
//	if((sim_params->tmask >> 15) & 0x1)fprintf(sim_params->outfile,"log(X): %f \t log(Evidence) estimate: %f\n",sim_params->logfactor*(sim_params->iter+1),sim_params->logZ); 
	if((sim_params->tmask >> 15) & 0x1)fprintf(sim_params->outfile,"log(X): %f \t log(Evidence) estimate: %f\n",sim_params->logX,sim_params->logZ); 
	if((sim_params->tmask >> 16) & 0x1)fprintf(sim_params->outfile,"Information estimate: %f\n",sim_params->H);
//	fprintf(stderr,"%f ",sim_params->logfactor * (sim_params->iter+1));
//	fprintf(stderr,"bad %f ",sim_params->Delta_logX * (sim_params->iter+1));
//	fprintf(stderr,",%f + %f * (%d - %d + 1) = %f, ",sim_params->logX_start,sim_params->Delta_logX,sim_params->iter,sim_params->iter_start,sim_params->logX);
	fprintf(stderr,"%f ",sim_params->logX);
#ifdef PARALLEL
	freemem_chain(&temporary);
    }
#endif
}



/* Nested Sampling routine.
   Eventually, we want to remove the reading in routine from here and put it into main. */
void nestedsampling(int thinning, int maxiter, simulation_params *sim_params){

  //set P and rank for serial job
  int P = 1; //number of processors used for NS
  int rank = 0; //master process
  void *comm_pointer = NULL;

  //NS points and their chaint-s
  Chain* cpoints = NULL;
  Chaint *chaint= (Chaint*)malloc(sizeof(Chaint));
  chaint->aat = NULL; chaint->ergt = NULL;  chaint->xaat = NULL; chaint->xaat_prev = NULL;
  //biasmap
  Biasmap *biasmap = (Biasmap*)malloc(sizeof(Biasmap));
  biasmap->distb = NULL;
  //temporary chain for reading in
  Chain* temporary = NULL;
  temporary = (Chain *)malloc(sizeof(Chain));
  temporary->aa = NULL; temporary->xaa = NULL; temporary->erg = NULL; temporary->xaa_prev = NULL;
  //number of chains
  int N = 0;
  //number of chains stored on this processor
  int current_stored = 0;

  time_t timer1 = 0;
  time_t timer2 = 0;
#ifdef PARALLEL
     timer1 = MPI_Wtime();
#else
     timer1 = time(NULL);
#endif


#ifdef PARALLEL

  MPI_Comm FLEX_WORLD, NS_WORLD;
  int S = sim_params->flex_params.number_of_processors;
  int SP,Srank, NSgrouprank,FLEXgrouprank;

  /*If only reading in from checkpoint file then do not use FLEX processor, even if it was asked for*/
  if(sim_params->num_NS_per_checkpoint == -1 && S != 0){
    S=0;
    fprintf(stderr,"Note, ignoring FLEX options as -C -1,FILENAME is used!\n");
  }


  int in_NS = 1;
  setup_communicators(&FLEX_WORLD,&NS_WORLD,S,&P,&rank,&SP,&Srank,&NSgrouprank,&FLEXgrouprank,&in_NS);
  /* P and rank now refer to the rank and number of processors in NS_WORLD and
   * SP and srank refer to the rank and number of processors in FLEX_WORLD*/

  if(rank == 0){
    fprintf(stderr,"Nested Sampling using %d processors + %d FLEX processors\n",P,S);
  }

  if(in_NS == -1){
    ns_for_flex_processor(FLEX_WORLD,FLEXgrouprank,biasmap,sim_params);
    return;
  }
  else{
    comm_pointer = (void*)(&NS_WORLD);
  }
  //parallel variable changes
  maxiter /= P;
  maxiter++;
  thinning /= P;
  if(thinning == 0) thinning = 1;

#endif
  Instructions *instructions=NULL;
  int copies = 0;
  //in case restarting from checkpoint
  int only_output_checkpoint = 0; 
  if(sim_params->num_NS_per_checkpoint == -1)
    only_output_checkpoint = 1;
  else sim_params->num_NS_per_checkpoint /= P;
  if(sim_params->num_NS_per_checkpoint == 0) sim_params->num_NS_per_checkpoint = 1;



  //initialising Nested Sampling variables
  sim_params->iter_start = 0; //this is to make sure we restart after checkpoint at correct place 
  sim_params->iter = 0;
  sim_params->H = 0.0; //information estimate
  sim_params->logZ = -DBL_MAX; //evidence or partition function estimate
  sim_params->amplitude = -M_PI; //MCMC crankshaft amplitude
  sim_params->alpha = 1; //shrinkage ratio of X, the available prior space ratio
  sim_params->logX = 0; //log of the available prior space ratio estimate
  sim_params->logX_start = 0; //log of the available prior space ratio estimate -- this is important if we restart the simulation
  sim_params->Delta_logX = 0; //step length along the logX axis: logX_{i-1} - logX_i = log alpha
  sim_params->log_DeltaX = 0; //current width of the prior space weight: log( X_{i-1} - X_i ) = i * log(1 - alpha)

  //local Nested Sampling variables
  double currE, logZnew, lweight;

  //the array where the chains are stored 
  //this stores the processor, index and ll of cpoints arrays from ALL processors
  ChainHash* chainhash = NULL; //only used by master processor of parallel
  //this is used to bring in and send out the P chains at each NS interation and for new sample point generation by MCMC
  Chain* chaincopies = NULL;



  //read in chains from either PDB file or checkpoint file
  if(sim_params->restart_from_checkpoint == 1){
    //note if only outputting from checkpoint then the tests are performed in here and no chains are stored or sent
    current_stored = read_in_from_checkpoint(sim_params, &biasmap, temporary, &cpoints, &chaint, &chainhash, rank, P,comm_pointer,only_output_checkpoint);
  } else {
    //this will also map and initialise the chains
    time_t timer_init1 = 0, timer_init2 = 0;
#ifdef PARALLEL
     timer_init1 = MPI_Wtime();
#else
     timer_init1 = time(NULL);
#endif
    current_stored = read_in_from_pdb(sim_params, &biasmap, temporary, &cpoints, &chaint, &chainhash, rank, P,comm_pointer);
#ifdef PARALLEL
     timer_init2 = MPI_Wtime();
#else
     timer_init2 = time(NULL);
#endif
    if (rank == 0) {
      fprintf(stderr,"NS reading in from checkpoint time: %g\n",(double)timer_init2-(double)timer_init1);
    }
    initialize_all_pdb_chains(sim_params, &biasmap, temporary, &cpoints, &chaint, &chainhash, rank, current_stored, P,comm_pointer);
#ifdef PARALLEL
     timer_init1 = MPI_Wtime();
#else
     timer_init1 = time(NULL);
#endif
    if (rank == 0) {
      fprintf(stderr,"NS pdb initialising time: %g\n",(double)timer_init1-(double)timer_init2);
    }
  }
  N = sim_params->N;
#ifdef PARALLEL
  //fprintf(stderr,"BEFORE rank %d sim_params->sequence: %s sim_params->seq: %s sim_params->NAA: %d sim_params->Nchains: %d\n",rank,sim_params->sequence,sim_params->seq, sim_params->NAA, sim_params->Nchains);
  //Copy sequence from rank 0 to all other nodes, this will be needed when setting up MC move lookup table
  MPI_Bcast(&sim_params->NAA,1,MPI_INT,0,NS_WORLD);
  MPI_Bcast(&sim_params->Nchains,1,MPI_INT,0,NS_WORLD);
  //fprintf(stderr,"BEFORE rank %d sim_params->sequence: %s sim_params->seq: %s sim_params->NAA: %d sim_params->Nchains: %d\n",rank,sim_params->sequence,sim_params->seq, sim_params->NAA, sim_params->Nchains);
  if (rank != 0) {
	sim_params->seq = (char *)realloc(sim_params->seq,(sim_params->NAA+1) * sizeof(char));
	sim_params->sequence = (char *)realloc(sim_params->sequence,(sim_params->NAA+sim_params->Nchains) * sizeof(char));
  }
  MPI_Bcast(sim_params->seq,sim_params->NAA+1,MPI_CHAR,0,NS_WORLD);
  MPI_Bcast(sim_params->sequence,sim_params->NAA+sim_params->Nchains,MPI_CHAR,0,NS_WORLD);
  //fprintf(stderr,"AFTER rank %d sim_params->sequence: %s sim_params->seq: %s\n",rank,sim_params->sequence,sim_params->seq);
#endif

  //for debugging, output checkpoint file
  //simulation_params temp_sim_params;
  //sim_params_copy(&temp_sim_params, sim_params);
  ////temp_sim_params.iter = 0;
  //temp_sim_params.num_NS_per_checkpoint = 1;
  //temp_sim_params.checkpoint_filename = "temp_check";
  //output_checkpoint_file(cpoints, current_stored, &temp_sim_params,comm_pointer);
  //stop("OK");

  //if we do not want to run a NS quit
  if(only_output_checkpoint == 1){
    if (rank == 0) {
      //do the tests and then free up and return  //for (int i=0; i<current_stored; i++) {      //  tests(&cpoints[i],biasmap,sim_params->tmask, sim_params, 0x11);      //}
      //fprintf(stderr,"ERROR ERROR ERROR, only outputting the chains stored on processor 0, bug needs to be fixed\n"); fflush(stderr);
      freemem_chaint(chaint);
      free(chaint);
      free(chainhash);
      biasmap_finalise(biasmap);
      freemem_chain(temporary);
    }
    //for (int i=0; i<current_stored; i++) {   //  freemem_chain(&cpoints[i]);   //}
    free(cpoints);
    return;
  }

  //set the checkpoint counter and initialise the output file
  if (rank == 0) {
    //The checkpoint counter needs advancing not to overwrite the checkpoint file that's just been read in
    if(sim_params->restart_from_checkpoint == 1){
      sim_params->checkpoint_counter++;
    }

    //initialise the new output file
    char *out = (char*)malloc(sizeof(char)*1010);
    if(sim_params->outfile_name != NULL){
      sprintf(out,"%s_%d",sim_params->outfile_name,sim_params->checkpoint_counter);
      if (sim_params->outfile) fclose(sim_params->outfile);
      sim_params->outfile = fopen(out, "w");
    }
    free(out);
  }

#ifdef PARALLEL
  //receive final values from rank 0
  MPI_Bcast(&N,1,MPI_INT,0,NS_WORLD);
  MPI_Bcast(&sim_params->iter_start,1,MPI_INT,0,NS_WORLD);
  MPI_Bcast(&(sim_params->checkpoint_counter),1,MPI_INT,0,NS_WORLD);
  MPI_Bcast(&(sim_params->amplitude),1,MPI_DOUBLE,0,NS_WORLD);
  if(sim_params->number_initial_MC != 0)  MC_first(chainhash,cpoints, chaint,  current_stored, biasmap,sim_params,  rank,  N, P, &NS_WORLD);
#endif


  /* register sample points with their logL in the heap */
  if (rank == 0) {
    constructhashheap(chainhash,N);
  }


  //setup chaincopies for the chains that will be collected every NS iteration to generate the new sample points, P for parallel, 1(=P) for serial
  //the slave processors will have 1 chaincopy, with the one they do MCMC on
  int size_of_chaincopies = 1;
  int size_of_instruction_set = 0;
  instructions = NULL;
#ifndef FAST
  if (rank == 0) {
    size_of_chaincopies = P;
  }
#endif
  chaincopies = (Chain*)malloc(sizeof(Chain)*size_of_chaincopies);
  for(int i = 0; i < size_of_chaincopies; i++){
    chaincopies[i].aa = NULL; chaincopies[i].erg = NULL;
    chaincopies[i].NAA = temporary->NAA; chaincopies[i].Nchains = temporary->Nchains;
    chaincopies[i].xaa = NULL; chaincopies[i].xaa_prev = NULL;
    allocmem_chain(&chaincopies[i],chaincopies[i].NAA,chaincopies[i].Nchains);
    for(int k = 0; k < temporary->NAA; k++){
      chaincopies[i].aa[k].id = temporary->aa[k].id;
      chaincopies[i].aa[k].num = temporary->aa[k].num;
      chaincopies[i].aa[k].etc = temporary->aa[k].etc;
    }
  }

#ifdef FAST
  size_of_instruction_set = 1;
  if(rank == 0){
    size_of_instruction_set = P;
  }
  instructions = (Instructions*)malloc(sizeof(Instructions)*size_of_instruction_set);
  for(int k = 0; k < size_of_instruction_set; k++ ){
    initialize_instruction_set(&(instructions[k]),P);
  }
#endif

  //End of reading and storing chains
  //Set up nested sampling


  //NS variables
  //shrinakge ratio
  if (rank == 0) {
#ifdef PARALLEL
    sim_params->alpha = 1.0 - (double)(P)/(double)(N+1);
#else
    //do not approximate exp(-1/N) as N/(N+1)
    sim_params->alpha = exp(-1.0/(double)N);
#endif

    //Delta logX
    sim_params->Delta_logX = log( sim_params->alpha);
    //log DeltaX
    sim_params->log_DeltaX = log ( 1.0 - sim_params->alpha); //log(P) - log(N+1); for parallel

    //adjust log_DeltaX (the width of the prior space) if we are restarting
    if(sim_params->restart_from_checkpoint == 1){
      sim_params->log_DeltaX += sim_params->logX_start;
    }
  }


  int converged = 1;

#ifdef PARALLEL
  if(rank == 0 && S != 0){
    /*send generic chain details to FLEX_WORLD */
    MPI_Bcast(&chaincopies[0].NAA,1,MPI_INT,0,FLEX_WORLD);
    int i;
    for(i = 1; i < chaincopies[0].NAA; i++){
      MPI_Bcast(&chaincopies[0].aa[i].id,1,MPI_CHAR,0,FLEX_WORLD);
      MPI_Bcast(&chaincopies[0].aa[i].num,1,MPI_INT,0,FLEX_WORLD);
      MPI_Bcast(&chaincopies[0].aa[i].etc,1,MPI_INT,0,FLEX_WORLD);
    }
  }
#endif

#ifdef PARALLEL
     timer2 = MPI_Wtime();
#else
     timer2 = time(NULL);
#endif
  if (rank == 0) {
     fprintf(stderr,"NS start-up time: %g\n",(double)timer2-(double)timer1);
  }


  /* MAIN NESTED SAMPLING LOOP */

  for(sim_params->iter = sim_params->iter_start; sim_params->iter < maxiter && converged; sim_params->iter++){
    check_to_output_checkpoint_file( cpoints,current_stored,N,sim_params,  rank, comm_pointer);

    int heaplength = N; //only master uses it

    /*if(rank == 0){
      for(int i = 1; i <=N; i++){
        fprintf(stderr,"CH %d: %d %d %lf\n",i,chainhash[i].processor,chainhash[i].index, chainhash[i].ll);
      }
    }

    for(int i = 0; i < current_stored; i++){
      fprintf(stderr,"CP: rank %d  index %d %lf %lf\n",rank,i,cpoints[i].ll,-totenergy(&(cpoints[i])));
    }*/



    //find the worst P samples (1 in case of serial)
    //and set logLstar accordingly
    //update all NS parameters
    if (rank == 0) {
      find_worst(sim_params,chainhash, &heaplength,  N, P);
      update_NS_parameters(sim_params, &converged, &logZnew,&lweight);
    }

    /*output sample point */
    //this works for both parallel (master/slave) and serial
    //chainhash and heap are set to NULL at the beginning for those who don't need it
    if((sim_params->iter+1) % thinning == 0){
      output_NS_point(chainhash, cpoints, biasmap, sim_params, N,comm_pointer);
    }

#ifdef PARALLEL
    collect_chains(chainhash,cpoints,chaincopies,sim_params,rank,P,N,&NS_WORLD, instructions);
#else
    do
      copies = 1 + ((int)(rand()/(double)RAND_MAX * N)) % N;
    while(copies == 1 && N > 1);
    copybetween(&chaincopies[0],&cpoints[chainhash[copies].index]);
#endif


    //explore constrained space, set logLstar
    currE = -chaincopies[0].ll;       

    //generate new points
    //stepsize();
    for(int i = 0; i < sim_params->iter_max; i++){
      move(&(chaincopies[0]),chaint,biasmap,sim_params->logLstar,&currE,0, sim_params);
    }
    //stepsize(); fprintf(sim_params->outfile,"Energy = %f Log(X) = %f\n",-logLstar,(sim_params->iter+1)*sim_params->Delta_logX);

    //update logLikelihood
    chaincopies[0].ll = -currE;



    //collect new points on the main processor in the heap
#ifdef PARALLEL
    return_and_reheap_chains(chainhash,cpoints,chaincopies,sim_params,rank,P,&heaplength,&NS_WORLD, instructions);
#else
    //put the new one back in heap and send them to overwrite the worst ones
    for(int k = 0; k < P; k++){
      heaplength++;
      chainhash[heaplength].ll = chaincopies[k].ll;
      copybetween(&cpoints[chainhash[heaplength].index],&chaincopies[k]);
      heapifyhashin(chainhash,heaplength);
    }

#endif



    //shrink interval
    if (rank == 0) {
      sim_params->log_DeltaX += sim_params->Delta_logX;
    }


    //time to change amplitude
    //use every chain for parallel simulations, and 5 for serial runs
#ifdef PARALLEL
    if((sim_params->iter*2*P)%N == 0){
#else
      if(sim_params->iter % N == 0){
#endif
        new_amplitude(&cpoints, biasmap, current_stored, sim_params,comm_pointer);
      }


      //has Z converged and deal with flex if required
#ifdef PARALLEL
      MPI_Bcast(&converged,1,MPI_INT,0,NS_WORLD);
      if(S != 0){
        int finished = converged;
        if(sim_params->iter == maxiter - 1) finished = 0;
        deal_with_flex(N,finished,rank, sim_params, NS_WORLD, FLEX_WORLD,chainhash,temporary,cpoints);
      }
#endif



    } //end of main loop

#ifdef PARALLEL
     timer1 = MPI_Wtime();
#else
     timer1 = time(NULL);
#endif
  if (rank == 0) {
     fprintf(stderr,"NS loop time: %g\n",(double)timer1-(double)timer2);
  }

    //clean up
#ifdef FAST
    for(int k = 0; k < size_of_instruction_set; k++){
      finalize_instruction_set(&(instructions[k]));
    }
    free(instructions);
#endif
    if (rank == 0) {
      //clean up rank 0 only memory, see later for rest of cleanup
      free(chainhash);
    }
    for(int i = 0; i < size_of_chaincopies; i++){
      freemem_chain(&(chaincopies[i]));
    }
    free(chaincopies);

    freemem_chain(temporary); free(temporary);
    freemem_chaint(chaint);
    free(chaint);
    biasmap_finalise(biasmap);
    for(int i = 0; i < current_stored; i++){
      freemem_chain(&cpoints[i]);
    }
    free(cpoints);

    return;


  }
