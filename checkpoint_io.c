/*
** Checkpoint input-output routines.
**
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>
#include<string.h>
#include<time.h>
#ifdef PARALLEL
#include<mpi.h>
#include"random16.h"
#endif
#include"error.h"
#include"params.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"
#include"vdw.h"
#include"energy.h"
#include"metropolis.h"
#include"probe.h" 
#include"checkpoint_io.h"

//====================================================================//
//                                                                    //
//                      PRINTING CHECKPOINT FILE                      //
//                                                                    //
//====================================================================//

/* Initialise the checkpoint file as sim_params->checkpoint_file,
   to be read in (action="r") or written (action="w" or "a"). */
void open_next_checkpoint_file(simulation_params *sim_params,char *action){

  char *checkpoint_name = (char*)malloc(sizeof(char)*1010);

  if(!sim_params->outfile_name)
     stop("Cannot open checkpoint file. Output name was not given.");

  //put together the name of the checkpoint file
  sprintf(checkpoint_name,"%s_%d",sim_params->checkpoint_filename,sim_params->checkpoint_counter);
  fprintf(stderr, "Checkpoint file: %s\n",checkpoint_name);

  //try to open it
  sim_params->checkpoint_file = fopen(checkpoint_name, action);
  if(!sim_params->checkpoint_file){
     stop("Error, checkpoint file cannot be opened");
  }

//  fprintf(stderr,"Successfully initialized checkpoint file\n");	
//  fprintf(stderr,"checkpoint counter %d\n",sim_params->checkpoint_counter);
  free(checkpoint_name);

}


/* Print header into the current checkpoint file.
   The header contains information about the peptide chains (NAA, N, seq)
   and NS simulation variables (iter_start, logL*, logZ, H, amplitude) */
//!! current_iter (sim_params->iter) is how many iterations the master processor has done,
//!! for serial runs it will be num_proc times larger than for parallel runs.
void print_checkpoint_header(simulation_params *sim_params, FILE *outfile, int current_iter){

  // print numbers of aa-s and NS points
  fprintf(outfile,"%d %d %d\n",sim_params->NAA,sim_params->N,sim_params->Nchains);
  // print and set sequence
  fprintf(outfile,"%s\n",sim_params->seq);

  // print global parameters and counters
  fprintf(outfile,"%d %lf %lf %lf %lf %lf\n", current_iter,
				sim_params->logX,
				sim_params->logLstar,
				sim_params->logZ,
				sim_params->H,
				sim_params->amplitude);

}


/* Print a checkpoint entry, that is a peptide chain. */
void print_checkpoint_entry(Chain *cpoints, simulation_params *sim_params, FILE *outfile, int N){
  int chainloop, aaloop,i;
  
  //fprintf(stderr,"Printing N=%d molecules\n",N);
  for(chainloop = 0; chainloop < N; chainloop++){
	 aaloop = 0;
	/* print first xaa (CA-CA) vector */
	 for(i = 0; i < 3; i++){
	    fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].xaa[aaloop][i][0], cpoints[chainloop].xaa[aaloop][i][1], cpoints[chainloop].xaa[aaloop][i][2]);
      }
	/* loop over the amino acids */	
	for(aaloop = 1; aaloop < cpoints[chainloop].NAA; aaloop++){
	/* print the atomic coordinates */
	  if(cpoints[chainloop].aa[aaloop].id != 'P'){
	    fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].h[0],cpoints[chainloop].aa[aaloop].h[1],cpoints[chainloop].aa[aaloop].h[2]);
	  }
	  fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].n[0],cpoints[chainloop].aa[aaloop].n[1],cpoints[chainloop].aa[aaloop].n[2]);
	  fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].ca[0],cpoints[chainloop].aa[aaloop].ca[1],cpoints[chainloop].aa[aaloop].ca[2]);
	  fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].c[0],cpoints[chainloop].aa[aaloop].c[1],cpoints[chainloop].aa[aaloop].c[2]);
	  fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].o[0],cpoints[chainloop].aa[aaloop].o[1],cpoints[chainloop].aa[aaloop].o[2]);
	  if(cpoints[chainloop].aa[aaloop].id != 'G'){
	    fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].cb[0],cpoints[chainloop].aa[aaloop].cb[1],cpoints[chainloop].aa[aaloop].cb[2]);
	    if(((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) && (cpoints[chainloop].aa[aaloop].id != 'A')){
	      fprintf(outfile,"%12.8f %12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].g[0],cpoints[chainloop].aa[aaloop].g[1],cpoints[chainloop].aa[aaloop].g[2],cpoints[chainloop].aa[aaloop].chi1);
	        if(cpoints[chainloop].aa[aaloop].id == 'V' || cpoints[chainloop].aa[aaloop].id == 'I' || cpoints[chainloop].aa[aaloop].id == 'T'){
	        fprintf(outfile,"%12.8f %12.8f %12.8f %12.8f\n", cpoints[chainloop].aa[aaloop].g2[0],cpoints[chainloop].aa[aaloop].g2[1],cpoints[chainloop].aa[aaloop].g2[2],cpoints[chainloop].aa[aaloop].chi2);
	      }
	    }
	  }
	    /* print next xaa (CA-CA) vector */
	  for(i = 0; i < 3; i++){
	    fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].xaa[aaloop][i][0], cpoints[chainloop].xaa[aaloop][i][1], cpoints[chainloop].xaa[aaloop][i][2]);
      }

	    /* print etc and the number of amino acid in the sequence */
	  fprintf(outfile,"%x %d %d\n",cpoints[chainloop].aa[aaloop].etc,cpoints[chainloop].aa[aaloop].num,cpoints[chainloop].aa[aaloop].chainid);   
    }
    /* print xaa (CA-CA) vector for multiple chain starts */
    for (int chainid = 0; chainid <= cpoints[chainloop].aa[cpoints[chainloop].NAA-1].chainid; chainid++) {
	 for(i = 0; i < 3; i++){
	    fprintf(outfile,"%12.8f %12.8f %12.8f\n", cpoints[chainloop].xaa_prev[chainid][i][0], cpoints[chainloop].xaa_prev[chainid][i][1], cpoints[chainloop].xaa_prev[chainid][i][2]);
	    //fprintf(stderr,"check point%d aa%d %12.8f %12.8f %12.8f\n", chainloop, chainid, cpoints[chainloop].xaa_prev[chainid][i][0], cpoints[chainloop].xaa_prev[chainid][i][1], cpoints[chainloop].xaa_prev[chainid][i][2]);
      }
    }
  
	/* print energy matrix */
    for(aaloop = 0; aaloop < cpoints[chainloop].NAA*cpoints[chainloop].NAA; aaloop++){
	  fprintf(outfile,"%12.8f ",cpoints[chainloop].erg[aaloop]);	   
    } 			
  }
}


/* Output checkpoint file for serial or parallel simulations,
   all-inclusive: from opening the file to printing the last line 
   The master processor writes its chains into a new file, then
   all other processors append it with their chains one after the other. */
void output_checkpoint_file(Chain *cpoints, int current_stored, simulation_params *sim_params,void *mpi_comm){

	//nested.c: j must be saved as sim_params->iter by now
	//	    add if(j != iter_start)
	//	    sim_params->NAA must be the same as cpoints[0].NAA
	//	    sim_params->N must be the same as N
	//	    sim_params->seq must be the same as seq
	//	    also: logLstar, logZ, H
	//	  what is tag used for?
	//	  current_stored needs updating before printing for each processor!

	//fprintf(stderr,"num_NS_per_checkpoint %d\n",sim_params->num_NS_per_checkpoint);
	if( sim_params->iter % sim_params->num_NS_per_checkpoint != 0 ) return;

	int rank = 0;
	char action[2];

#ifdef PARALLEL
	MPI_Comm *MPI_COMM = mpi_comm;
	int P = 1;
	int tag = 0;
	MPI_Status status;
	MPI_Comm_size(*MPI_COMM, &P);
	MPI_Comm_rank(*MPI_COMM, &rank);
#endif

//	fprintf(stderr,"j(sim_params->iter)= %d and sim_params->iter_start= %d\n",sim_params->iter,sim_params->iter_start);

#ifdef PARALLEL
	if (rank != 0) {
		//wait until it's my turn to output my cpoints
		MPI_Recv(&tag,1,MPI_INT,rank-1,0,*MPI_COMM,&status);
	}
#endif

	//The master processor will create the file, the others will append to it.
	if (rank == 0) {
		sprintf(action,"w");
	} else {
		sprintf(action,"a");
	}

	open_next_checkpoint_file(sim_params,action);

	if (rank == 0) {
		print_checkpoint_header(sim_params, sim_params->checkpoint_file, sim_params->iter);
	}

	print_checkpoint_entry(cpoints, sim_params, sim_params->checkpoint_file, current_stored);

	fclose(sim_params->checkpoint_file);
	sim_params->checkpoint_file = NULL;

	sim_params->checkpoint_counter++;

#ifdef PARALLEL
	if (rank == 0 && P > 1) {
		//tell process 1 it's his turn
		MPI_Send(&tag,1,MPI_INT,1,0,*MPI_COMM);
		//wait to hear from processor P-1 that we're done
		MPI_Recv(&tag,1,MPI_INT,P-1,0,*MPI_COMM,&status);
	} else {
		//tell next processor it's their turn (or tell 0 we're done if rank = P-1)
		MPI_Send(&tag,1,MPI_INT,(rank+1) % P,0,*MPI_COMM);
	}
#endif

}



//====================================================================//
//                                                                    //
//                      READING  CHECKPOINT FILE                      //
//                                                                    //
//====================================================================//



/* Open the checkpoint file and read in its header.  Update sim_params.
   The header contains information about the peptide chains (NAA, N, seq)
   and NS simulation variables (iter_start, logL*, logZ, H, amplitude) */
void read_checkpoint_header(simulation_params *sim_params){

  if(NULL == sim_params->checkpoint_file){
	stop("read_checkpoint_header: Checkpoint file is not open yet.");
  }

  //read numbers of aa-s and NS points
  int k = 0;
  if ((k = fscanf(sim_params->checkpoint_file,"%d %d %d\n",&(sim_params->NAA),&(sim_params->N),&(sim_params->Nchains))) != 3) {
	stop("read_checkpoint_header: Could not read amino acid and chain numbers.\n");
  }
  // read in and set sequence
  sim_params->seq = (char*)malloc(sizeof(char)*(sim_params->NAA+2));
  if ((k = fscanf(sim_params->checkpoint_file,"%s\n",sim_params->seq)) != 1) {
	stop("read_checkpoint_header: Could not read sequence.\n");
  }

  // read in global parameters and counters
  if ((k = fscanf(sim_params->checkpoint_file,"%d %lf %lf %lf %lf %lf\n", &(sim_params->iter_start),
				&(sim_params->logX_start),
				&(sim_params->logLstar),
				&(sim_params->logZ),
				&(sim_params->H),
				&(sim_params->amplitude))) != 6) {
	stop("read_checkpoint_header: Could not read global parameters and counters.\n");
  }
  sim_params->logX = sim_params->logX_start;
  fprintf(stderr,"Checkpoint file opened\n iter_start:%d logX:%f L*:%lf logZ:%lf (0.0 means -DBL_MAX) H:%lf amplitude:%lf\n",sim_params->iter_start,sim_params->logX_start,sim_params->logLstar,sim_params->logZ,sim_params->H,sim_params->amplitude);

  if((sim_params->logZ) == 0.0) (sim_params->logZ) = -DBL_MAX;

}


/* Read a peptide chain entry from the checkpoint file.
   The atomic coordinates are given in rows, then the energy matrix. */
void read_checkpoint_entry(Chain *cpoints, simulation_params *sim_params){

    int  aaloop,i;

	aaloop = 0;

	int k = 0;
	/* read in first xaa (CA-CA) vector */
	for(i = 0; i < 3; i++){
	    if ((k = fscanf(sim_params->checkpoint_file,"%lf %lf %lf\n", &(cpoints->xaa[aaloop][i][0]), &(cpoints->xaa[aaloop][i][1]), &(cpoints->xaa[aaloop][i][2]))) != 3) {
		stop("read_checkpoint_entry: Could not read first amino acid xaa vector.\n");
	    }
	}

	/* set the aa[0] (this is not stored in the checkpoint file) */
	cpoints->aa[0].id = 'A';
	cpoints->aa[0].num = 0;
	cpoints->aa[0].chainid = 0;
	cpoints->aa[0].etc = 0;

	/* loop over the amino acids */	
	for(aaloop = 1; aaloop < sim_params->NAA; aaloop++){
	/* read in the atomic coordinates */
	    if(cpoints->aa[aaloop].id != 'P'){
		if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->aa[aaloop].h[0]),&(cpoints->aa[aaloop].h[1]),&(cpoints->aa[aaloop].h[2]))) != 3) {
			stop("read_checkpoint_entry: Could not read amide hydrogen coordinates.\n");
		}
	    }
	    if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->aa[aaloop].n[0]),&(cpoints->aa[aaloop].n[1]),&(cpoints->aa[aaloop].n[2]))) != 3) {
		stop("read_checkpoint_entry: Could not read amide nitrogen coordinates.\n");
	    }
	    if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->aa[aaloop].ca[0]),&(cpoints->aa[aaloop].ca[1]),&(cpoints->aa[aaloop].ca[2]))) != 3) {
		stop("read_checkpoint_entry: Could not read alpha carbon coordinates.\n");
	    }
	    if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->aa[aaloop].c[0]),&(cpoints->aa[aaloop].c[1]),&(cpoints->aa[aaloop].c[2]))) != 3) {
		stop("read_checkpoint_entry: Could not read carbonil carbon coordinates.\n");
	    }
	    if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->aa[aaloop].o[0]),&(cpoints->aa[aaloop].o[1]),&(cpoints->aa[aaloop].o[2]))) != 3) {
		stop("read_checkpoint_entry: Could not read oxygen coordinates.\n");
	    }
		//TODO: we could invalidate CG atoms by setting chi1 or chi2 to DBL_MAX
	    if(cpoints->aa[aaloop].id != 'G'){
		if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->aa[aaloop].cb[0]),&(cpoints->aa[aaloop].cb[1]),&(cpoints->aa[aaloop].cb[2]))) != 3) {
		    stop("read_checkpoint_entry: Could not read beta carbon coordinates.\n");
		}
	        if(((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) && (cpoints->aa[aaloop].id != 'A')){
		    if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf %lf\n", &(cpoints->aa[aaloop].g[0]),&(cpoints->aa[aaloop].g[1]),&(cpoints->aa[aaloop].g[2]),&(cpoints->aa[aaloop].chi1))) != 4) {
			stop("read_checkpoint_entry: Could not read gamma side chain coordinates.\n");
		    }
		    if(cpoints->aa[aaloop].id == 'V' || cpoints->aa[aaloop].id == 'I' || cpoints->aa[aaloop].id == 'T'){
			if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf %lf\n", &(cpoints->aa[aaloop].g2[0]),&(cpoints->aa[aaloop].g2[1]),&(cpoints->aa[aaloop].g2[2]),&(cpoints->aa[aaloop].chi2))) != 4) {
			    stop("read_checkpoint_entry: Could not read V/I/T gamma2 side chain coordinates.\n");
			}
		    }
		}
	    }
	    /* read in next xaa (ca-ca) vector */
	    for(i = 0; i < 3; i++){
		if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->xaa[aaloop][i][0]), &(cpoints->xaa[aaloop][i][1]), &(cpoints->xaa[aaloop][i][2]))) != 3) {
		    stop("read_checkpoint_entry: Could not read xaa vector.\n");
		}
	    }
	    /* read in etc and the number of amino acid in the sequence */
	    if ((k = fscanf(sim_params->checkpoint_file, "%x %d %d\n",&(cpoints->aa[aaloop].etc),&(cpoints->aa[aaloop].num),&(cpoints->aa[aaloop].chainid))) != 3) {
		stop("read_checkpoint_entry: Could not read etc and NAA and chainid.\n");
	    }
	}
	/* read in the xaa (CA-CA) vector for multiple chain starts */
	for (int chainid = 0; chainid <= cpoints->aa[sim_params->NAA-1].chainid; chainid++) {
	    for(i = 0; i < 3; i++){
		if ((k = fscanf(sim_params->checkpoint_file, "%lf %lf %lf\n", &(cpoints->xaa_prev[chainid][i][0]), &(cpoints->xaa_prev[chainid][i][1]), &(cpoints->xaa_prev[chainid][i][2]))) != 3) {
		    stop("read_checkpoint_entry: Could not read xaa vector for multiple chain starts.\n");
		}
	    }
	}

	/* read in energy matrix */
	for(aaloop = 0; aaloop < sim_params->NAA * sim_params->NAA; aaloop++){
	    if ((k = fscanf(sim_params->checkpoint_file, "%lf ",&(cpoints->erg[aaloop]))) != 1) {
		stop("read_checkpoint_entry: Could not read energy matrix.\n");
	    }
	}

}


/* Read in all peptide chains from a checkpoint file, and distribute them into memory if running parallel. */
int read_in_from_checkpoint(simulation_params *sim_params,
			Biasmap **biasmap,
			Chain *temporary,
			Chain **cpoints,
			Chaint **chaint,
			ChainHash **chainhash,
			int rank,
			int P, void*mpi_comm, int only_output_checkpoint){
#ifdef PARALLEL
  MPI_Comm *MPI_COMM = mpi_comm;
#endif


    int i;
    int current_stored = 0;

    /* Set up the memory buffer on all nodes */
    /* peptide chains */
    *cpoints = NULL;
    /* peptide chain for MC moves */
    if (*chaint) freemem_chaint(*chaint);
    *chaint = (Chaint*)realloc(*chaint,sizeof(Chaint));
    (*chaint)->aat = NULL; (*chaint)->ergt = NULL;  (*chaint)->xaat = NULL; (*chaint)->xaat_prev = NULL;
    /* bias map */
    *biasmap = (Biasmap*)realloc(*biasmap,sizeof(Biasmap));
    (*biasmap)->distb = NULL;

    //temporary chain for reading
    temporary->aa = NULL; temporary->xaa = NULL; temporary->erg = NULL; temporary->xaa_prev = NULL;

    // only read on the master processor
    if(rank == 0){

	//this stores the processor, index and ll of cpoints arrays from ALL processors
	*chainhash = (ChainHash*)malloc(sizeof(ChainHash));
   
	int counter = 0; 

	//read general and global parameters from the checkpoint header
	open_next_checkpoint_file(sim_params,"r");
	read_checkpoint_header(sim_params);
//	N = sim_params->N;
	//fprintf(stderr,"checkpoint counter %d\n",sim_params->checkpoint_counter);
	//fprintf(stderr,"sequence %s\n",sim_params->seq);

	//initialise the temporary chain into which we read in all the configs
	temporary->NAA = sim_params->NAA;
	temporary->Nchains = sim_params->Nchains;
	allocmem_chain(temporary,sim_params->NAA,sim_params->Nchains);
	for (i=0; i< sim_params->NAA; i++) {
	    temporary->aa[i].id = sim_params->seq[i]; 
	    temporary->aa[i].num = i; 
	}
	//set aa[0] chain id, as it will not be read in.
	temporary->aa[0].chainid = 0;
	// init biasmap
	biasmap_initialise(temporary,*biasmap,&(sim_params->protein_model));
	aat_init(temporary,*chaint);


	fprintf(stderr,"reading in checkpoint snapshots...\n");

	int breakk = 0;
	//read in chains into the temporary chain, then save them
	while (breakk==0) {
	    //PDB file
	    if (sim_params->NAA != temporary->NAA) stop("nestedsampling: sim_params->NAA != temporary->NAA");
	    /* read in coordinates and energy matrix */
	    read_checkpoint_entry(temporary, sim_params);
	  
	    /* calc total energy */
	    temporary->ll = -totenergy(temporary); 

	    /* store details in main heap (indexed from 1 to be awkward) */
	    if(only_output_checkpoint != 1){
	      counter = store_chain(chainhash,temporary, *biasmap, *chaint, P,&current_stored, counter, cpoints, sim_params, rank,mpi_comm);
	    }
	    else{
	      tests(temporary,*biasmap,sim_params->tmask, sim_params, 0x11, NULL);
	      counter++;
	    }
	    if (counter >= sim_params->N) breakk=1;
	}

	fprintf(stderr,"finished reading.\n");
	//N = counter;
	/* update sim_params in case N changed */
        sim_params->N = counter;
	/* sequence contatining separators for multi-chain proteins */
	fprintf(stderr,"Saving sequence for multi-chain protein.\n");
	sim_params->sequence = (char *)realloc(sim_params->sequence,(sim_params->NAA+temporary->aa[sim_params->NAA-1].chainid) * sizeof(char));
	sim_params->sequence[0] = 'A';
	int next = 1;
	int Nchains = 1;
	for (int i=1; i< sim_params->NAA; i++) {
	    if (i > 1 && temporary->aa[i].chainid != temporary->aa[i-1].chainid) {
		sim_params->sequence[next] = '_';
		next ++;
		Nchains ++;
	    }
	    sim_params->sequence[next] = temporary->aa[i].id;
	    next ++;
	}
	sim_params->sequence[sim_params->NAA+temporary->aa[sim_params->NAA-1].chainid-1] = '\0';
	if (Nchains != temporary->aa[sim_params->NAA-1].chainid) {
		fprintf(stderr,"Nchains = %d, last_chainid = %d\n",Nchains,temporary->aa[sim_params->NAA-1].chainid);
		stop("read_in_from_checkpoint: The number of chains != the last chain ID.\n");
	}
	if (Nchains != temporary->Nchains) {
		fprintf(stderr,"Nchains = %d, nchains = %d\n",Nchains,temporary->Nchains);
		stop("read_in_from_checkpoint: The number of chains != number of chains in the read PDB.\n");
	}
	sim_params->Nchains = Nchains;

    }
#ifdef PARALLEL
    else {
	// receive and store all chains that processor 0 sends
      if(only_output_checkpoint != 1){
    	store_chain(NULL,temporary, *biasmap, *chaint, P, &current_stored, 0, cpoints, sim_params, rank,mpi_comm);
      }
    }

    if (rank == 0) {
	//tell other processors there are no more chains coming
	for(int i = 1; i < P; i++){
	    int minus1 = -1;
	    MPI_Send(&minus1,1,MPI_INT,i,current_stored,*MPI_COMM);
	}
    }
#endif

    if (rank==0) {
	//close checkpoint file (successfully read)   
	fclose(sim_params->checkpoint_file);
	sim_params->checkpoint_file = NULL;
    }
    if (rank==0) fprintf(stderr,"Successfully initialized checkpoint file\n");	

    return current_stored;
//all have been read in by now
//close checkpoint file in nestedsampling

}


/* Store chain in the memory and advance the counter by 1.
   For a serial program, all chains are stored on the only processor.
   For a parallel program, the chains are distributed across the processors,
   and the master processor keeps track of which processor has got which chain
   using chainhash. */
int store_chain(ChainHash **chainhash, Chain *temporary, Biasmap *biasmap, Chaint *chaint, int P, int *current_stored, int counter, Chain **cpoints, simulation_params *sim_params, int rank, void*mpi_comm){

#ifdef PARALLEL
  MPI_Comm *MPI_COMM = mpi_comm;
  //if this is first one then give heads up to other processors
	if((rank == 0 && counter == 0) || rank != 0){ //counter = 0 only holds for master
	    //send and get NAA, seq etc from processor 0
	    MPI_Bcast(&(temporary->NAA),1,MPI_INT,0,*MPI_COMM);
	    MPI_Bcast(&(temporary->Nchains),1,MPI_INT,0,*MPI_COMM);
	    if (rank != 0) {
		allocmem_chain(temporary,temporary->NAA,temporary->Nchains);
	    }
	    for(int i = 0; i < temporary->NAA; i++){
	        MPI_Bcast(&(temporary->aa[i].num),1,MPI_INT, 0,*MPI_COMM);
	        MPI_Bcast(&(temporary->aa[i].etc),1,MPI_INT, 0,*MPI_COMM);
	        MPI_Bcast(&(temporary->aa[i].id) ,1,MPI_CHAR,0,*MPI_COMM);
	        MPI_Bcast(&(temporary->aa[i].chainid),1,MPI_INT, 0,*MPI_COMM);
	    }  
	    if (rank != 0) {
		//initiate	
		biasmap_initialise(temporary,biasmap,&(sim_params->protein_model));
		aat_init(temporary, chaint);
	    }
	}

    if (rank == 0) {
	//which processor gets it?
	//store details in main heap (indexed from 1 to be awkward)
        *chainhash = (ChainHash*)realloc(*chainhash,(counter+2)* sizeof(ChainHash));
        (*chainhash)[counter+1].processor = counter % P;
        (*chainhash)[counter+1].ll = temporary->ll; 
#else
	//store details in main heap (indexed from 1 to be awkward)
        *chainhash = (ChainHash*)realloc(*chainhash,(counter+2)* sizeof(ChainHash));
        (*chainhash)[counter+1].processor = counter % P; //=0 always on the master processor
        (*chainhash)[counter+1].ll = temporary->ll; 
#endif

#ifdef PARALLEL
        if(counter % P == 0) {
#endif
	    (*current_stored)++;
#ifdef PARALLEL
	}
        (*chainhash)[counter+1].index = *current_stored-1;
#else
        (*chainhash)[counter+1].index = *current_stored-1;
#endif
#ifdef PARALLEL

	//if keeping chain this processor...
	//this is the same as in serial
	if(counter % P == 0){
#endif
	    /* store on the master processor */
	    counter++;
	    *cpoints = (Chain*)realloc(*cpoints,*current_stored * sizeof(Chain));
	    (*cpoints)[*current_stored-1].NAA = temporary->NAA;
	    (*cpoints)[*current_stored-1].Nchains = temporary->Nchains;
	    (*cpoints)[*current_stored-1].aa = NULL;
	    (*cpoints)[*current_stored-1].xaa = NULL;
	    (*cpoints)[*current_stored-1].erg = NULL;
	    (*cpoints)[*current_stored-1].xaa_prev = NULL;
	    allocmem_chain(&((*cpoints)[*current_stored-1]),temporary->NAA,temporary->Nchains);
	    copybetween(&((*cpoints)[*current_stored-1]),temporary);
		//or sending to be stored on another processor...
#ifdef PARALLEL
	 } else{
		  int proc_to_send_to = counter % P;
		  counter++;
		  //tell processor it will be sent a chain
		  MPI_Send(current_stored,1,MPI_INT,proc_to_send_to,*current_stored,*MPI_COMM);
		  mpi_send_chain(temporary, 0, proc_to_send_to, &(sim_params->logLstar) , *current_stored,*MPI_COMM);
	}
    } else { // slave processor receiving chain
    
	MPI_Status status;

	//keep accepting new chains to store in cpoints until get sent -1
	int got_all = 0;
	while(got_all != -1){
	    MPI_Recv(&got_all,1,MPI_INT,0,MPI_ANY_TAG,*MPI_COMM,&status);
	    if(got_all != -1){ 
		(*current_stored)++;
		*cpoints = (Chain*)realloc(*cpoints,*current_stored * sizeof(Chain));
		(*cpoints)[*current_stored-1].NAA = temporary->NAA;
		(*cpoints)[*current_stored-1].Nchains = temporary->Nchains;
		(*cpoints)[*current_stored-1].aa = NULL;
		(*cpoints)[*current_stored-1].xaa = NULL;
		(*cpoints)[*current_stored-1].xaa_prev = NULL;
		(*cpoints)[*current_stored-1].erg = NULL;
		allocmem_chain(&((*cpoints)[*current_stored-1]),temporary->NAA,temporary->Nchains);
		mpi_rec_chain(&((*cpoints)[*current_stored-1]), 0, rank, &(sim_params->logLstar) , *current_stored,*MPI_COMM);
	    }
	}
    
    }
#endif
    return counter;

}



//====================================================================//
//                                                                    //
//                          READING PDB FILE                          //
//                                                                    //
//====================================================================//



/* Read in all peptide chains from a PDB file, and distribute them into memory if running parallel. */
int read_in_from_pdb(simulation_params *sim_params,
			Biasmap **biasmap,
			Chain *temporary,
			Chain **cpoints,
			Chaint **chaint,
			ChainHash **chainhash,
			int rank,
			int P, void *mpi_comm){

    int current_stored = 0;

    if (sim_params->restart_from_checkpoint) stop("Attempted to read in pdb when restarting from checkpoint file was requested.");

#ifdef PARALLEL
  MPI_Comm *MPI_COMM = mpi_comm;
#endif

//    if(rank == 0){
//	sim_params->infile = fopen(sim_params->infile_name,"r");
//    }

    /* Set up the memory buffer on all nodes */
    /* peptide chains */
    *cpoints = NULL;
    /* peptide chain for MC moves */
    if (*chaint) freemem_chaint(*chaint);
    *chaint = (Chaint*)realloc(*chaint,sizeof(Chaint));
    (*chaint)->aat = NULL; (*chaint)->ergt = NULL;  (*chaint)->xaat = NULL; (*chaint)->xaat_prev = NULL;
    /* bias map */
    *biasmap = (Biasmap*)realloc(*biasmap,sizeof(Biasmap));
    (*biasmap)->distb = NULL;

    //temporary chain for reading
    temporary->aa = NULL; temporary->xaa = NULL; temporary->erg = NULL; temporary->xaa_prev = NULL;
    temporary->ll = 0;

    // only read on the master processor
    if(rank == 0){

	//this stores the processor, index and ll of cpoints arrays from ALL processors
	*chainhash = (ChainHash*)malloc(sizeof(ChainHash));

	int counter = 0; 

	fprintf(stderr,"reading in PDB snapshots...\n");

	int breakk = 0;
	double time_readpdb = 0, time_storechain = 0;
	/*now read in pdbs from pdbfile */
	//read in chains into the temporary chain, then save them
	while (breakk==0) {
	    //PDB file
	    //read in the PDB and set up the biasmap and energy matrix
	    time_t timer1, timer2;
#ifdef PARALLEL
     timer1 = MPI_Wtime();
#else
     timer1 = time(NULL);
#endif
	    //read next PDB snapshot
		fprintf(stderr,"reading next chain...\n");
	    if (pdbin(temporary,sim_params,sim_params->infile) == EOF) {
		breakk = 1;
		continue;
	    }
		fprintf(stderr,"next chain: NAA=%d Nchains=%d\n",temporary->NAA,temporary->Nchains);
#ifdef PARALLEL
     timer2 = MPI_Wtime();
#else
     timer2 = time(NULL);
#endif
	    if (rank == 0) {
	      time_readpdb += (double)timer2 - (double)timer1;
	    }
	    //read biasmap and set sequence when reading the first snapshot
	    if(counter == 0){
		//initialize energy
		biasmap_initialise(temporary,*biasmap,&(sim_params->protein_model));
		char * seq = (char*)malloc(sizeof(char)*(temporary->NAA+1));
		seq[0] = 'A';
		for(int i = 1; i < temporary->NAA; i++) seq[i] = temporary->aa[i].id; 
		seq[temporary->NAA] = '\0'; 	 
		if (sim_params->seq) free(sim_params->seq);
		copy_string(&(sim_params->seq),seq);
		free(seq);
	    }
	  
	    // save chain
	    counter = store_chain(chainhash,temporary, *biasmap, *chaint, P, &current_stored, counter, cpoints, sim_params, rank,mpi_comm);
#ifdef PARALLEL
     timer1 = MPI_Wtime();
#else
     timer1 = time(NULL);
#endif
	    if (rank == 0) {
	      time_storechain += (double)timer1 - (double)timer2;
	    }

	}

	fprintf(stderr,"finished reading.\n");
	fprintf(stderr,"time taken for pdb reading: %g.\n",time_readpdb);
	fprintf(stderr,"time taken for pdb storing: %g.\n",time_storechain);
	//TODO: by now, temporary may have been corrupted in the last pdbin step, so use the first snapshot stored on rank 0

        sim_params->N = counter;
	sim_params->NAA = temporary->NAA;
	sim_params->seq = (char *)realloc(sim_params->seq,(sim_params->NAA+1) * sizeof(char));
	sim_params->seq[0] = 'A';
	for (int i=1; i< sim_params->NAA; i++) {
	    sim_params->seq[i] = temporary->aa[i].id; 
	}
	sim_params->seq[sim_params->NAA] = '\0';
	/* sequence contatining separators for multi-chain proteins */
	sim_params->sequence = (char *)realloc(sim_params->sequence,(sim_params->NAA+temporary->aa[sim_params->NAA-1].chainid) * sizeof(char));
	sim_params->sequence[0] = 'A';
	int next = 1;
	int Nchains = 1;
	for (int i=1; i< sim_params->NAA; i++) {
	    if (i > 1 && temporary->aa[i].chainid != temporary->aa[i-1].chainid) {
		sim_params->sequence[next] = '_';
		next ++;
		Nchains ++;
	    }
	    sim_params->sequence[next] = temporary->aa[i].id;
	    next ++;
	}
	sim_params->sequence[sim_params->NAA+temporary->aa[sim_params->NAA-1].chainid-1] = '\0';
	if (Nchains != temporary->aa[sim_params->NAA-1].chainid) {
		fprintf(stderr,"Nchains = %d, last_chainid = %d\n",Nchains,temporary->aa[sim_params->NAA-1].chainid);
		stop("read_in_from_pdb: The number of chains != the last chain ID.\n");
	}
	if (Nchains != temporary->Nchains) {
		fprintf(stderr,"Nchains = %d, nchains = %d\n",Nchains,temporary->Nchains);
		stop("read_in_from_pdb: The number of chains != number of chains in the read PDB.\n");
	}
	sim_params->Nchains = Nchains;

    }
#ifdef PARALLEL
    else { //rank != 0
	// receive and store all chains that processor 0 sends
	store_chain(NULL,temporary, *biasmap, *chaint, P, &current_stored, 0, cpoints, sim_params, rank,mpi_comm);
    }

    if (rank == 0) {
	//tell other processors there are no more chains coming
	for(int i = 1; i < P; i++){
   	    int minus1 = -1;
	    MPI_Send(&minus1,1,MPI_INT,i,current_stored,*MPI_COMM);
	}
    }
#endif

//    if(rank == 0){
//	fclose(sim_params->infile);
//	sim_params->infile = NULL;
//    }

    return current_stored;

}


/* Map all peptide chains onto the crankite model and calculate the energy matrices.
   Do this on all nodes if running parallel. */
void initialize_all_pdb_chains(simulation_params *sim_params,
			Biasmap **biasmap,
			Chain *temporary,
			Chain **cpoints,
			Chaint **chaint,
			ChainHash **chainhash,
			int rank,
			int current_stored,
			int P, void *mpi_comm){

#ifdef PARALLEL
    MPI_Comm *MPI_COMM = mpi_comm;
    MPI_Status status;
#endif

    //initialising all chains stored on all processors
    for (int i=0; i<current_stored; i++) {

	// fix peptide if needed
	if(sim_params->protein_model.fixit) fixpeptide((*cpoints)[i].aa, (*cpoints)[i].NAA, &(sim_params->protein_model));

	//mark constrained and fixed amino acids
	mark_fixed_aa_from_file(&((*cpoints)[i]),sim_params);
	mark_constrained_aa_from_file(&((*cpoints)[i]),sim_params);

	chkpeptide((*cpoints)[i].aa, (*cpoints)[i].NAA, &(sim_params->protein_model));
	//map chain onto the CRANKITE model
	initialize(&((*cpoints)[i]),*chaint,sim_params);

	//calculate the energy matrix
	energy_matrix_calculate(&((*cpoints)[i]),*biasmap,&(sim_params->protein_model));
	(*cpoints)[i].ll = -totenergy(&((*cpoints)[i]));

    }


    //collect all total energies for the ChainHash on the master processor
    if (rank==0) {
	for (int i=1; i<=sim_params->N; i++) {
	    int from=(*chainhash)[i].processor;
    	    if (from==0) { //copy from the master
		(*chainhash)[i].ll = (*cpoints)[(*chainhash)[i].index].ll;
	    } else { //copy from a slave
#ifdef PARALLEL
		int tag = (*chainhash)[i].index;
		MPI_Recv(&((*chainhash)[i].ll),1,MPI_DOUBLE,from,tag,*MPI_COMM,&status);
#else
		stop("In serial version the chains cannot be stored on processors other than the main one!");
#endif
	    }
	}
#ifdef PARALLEL
    //send total energies to the chainhash on the main processor
    } else {
	for (int i=0; i<current_stored; i++) {
	    int tag = i;
	    MPI_Send(&((*cpoints)[i].ll),1,MPI_DOUBLE,0,tag,*MPI_COMM);
	}
#endif
    }

}



//====================================================================//
//                                                                    //
//                         MPI  COMMUNICATION                         //
//                                                                    //
//====================================================================//


#ifdef PARALLEL
/* Send a peptide chain across MPI.
   Works with multi-chain proteins. */
void mpi_send_chain(Chain* nsconformation, int from, int to, double *logLstar, int iter, MPI_Comm MPI_COMM){


  int NAA = nsconformation->NAA;
  int Nchains = nsconformation->Nchains;
  int i,j, counter = 0;
  double coords[35*NAA];
  int etcs[NAA];
  int nums[NAA];
  char ids[NAA];
  int chainids[NAA];
  double xaa_prev[(Nchains+1)*9];
  if(from==0)MPI_Send(logLstar,1,MPI_DOUBLE,to,iter*10,MPI_COMM);
  MPI_Send(nsconformation->erg,NAA*NAA,MPI_DOUBLE,to,iter*10+1,MPI_COMM);
  MPI_Send(&(nsconformation->ll),1,MPI_DOUBLE,to,iter*10+2,MPI_COMM);
  for(j = 0; j < NAA; j++){
    etcs[j] = nsconformation->aa[j].etc;
    nums[j] = nsconformation->aa[j].num;
    ids[j] = nsconformation->aa[j].id;
    chainids[j] = nsconformation->aa[j].chainid;
    for(i = 0; i < 3; i++){
      coords[counter++] = nsconformation->aa[j].h[i];		
	  coords[counter++] = nsconformation->aa[j].n[i];		
	  coords[counter++] = nsconformation->aa[j].ca[i];		
	  coords[counter++] = nsconformation->aa[j].c[i];		
	  coords[counter++] = nsconformation->aa[j].o[i];
	  coords[counter++] = nsconformation->aa[j].cb[i];
	  coords[counter++] = nsconformation->aa[j].g[i];
	  coords[counter++] = nsconformation->aa[j].g2[i];
	  coords[counter++] = nsconformation->xaa[j][i][0]; 
	  coords[counter++] = nsconformation->xaa[j][i][1]; 
	  coords[counter++] = nsconformation->xaa[j][i][2]; 
    }		
    coords[counter++] = nsconformation->aa[j].chi1;
    coords[counter++] = nsconformation->aa[j].chi2;
    
  }
  counter = 0;
  for(j = 0; j <= Nchains; j++){
    for(i = 0; i < 3; i++){
	  xaa_prev[counter++] = nsconformation->xaa_prev[j][i][0]; 
	  xaa_prev[counter++] = nsconformation->xaa_prev[j][i][1]; 
	  xaa_prev[counter++] = nsconformation->xaa_prev[j][i][2]; 
    }
  }
  MPI_Send(coords,NAA*35,MPI_DOUBLE,to,iter*10+3,MPI_COMM);
  MPI_Send(etcs,NAA,MPI_INT,to,iter*10+4,MPI_COMM);
  MPI_Send(nums,NAA,MPI_INT,to,iter*10+5,MPI_COMM);
  MPI_Send(ids,NAA,MPI_CHAR,to,iter*10+6,MPI_COMM);
  MPI_Send(chainids,NAA,MPI_INT,to,iter*10+7,MPI_COMM);
  MPI_Send(&Nchains,1,MPI_INT,to,iter*10+8,MPI_COMM);
  MPI_Send(xaa_prev,(Nchains+1)*9,MPI_DOUBLE,to,iter*10+9,MPI_COMM);
}

/* Receive a chain through MPI.
   Works with multi-chain proteins. */
void mpi_rec_chain(Chain *nsconformation, int from, int to, double *logLstar, int iter, MPI_Comm MPI_COMM){
  int NAA = nsconformation->NAA;
  int i,j, counter = 0;
  double coords[35*NAA];
  int etcs[NAA];
  int nums[NAA];
  char ids[NAA];
  int chainids[NAA];
  int Nchains;
  MPI_Status info;
  if(from==0) MPI_Recv(logLstar,1,MPI_DOUBLE,from,iter*10,MPI_COMM,&info);
  MPI_Recv(nsconformation->erg,NAA*NAA,MPI_DOUBLE,from,iter*10+1,MPI_COMM,&info);
  MPI_Recv(&(nsconformation->ll),1,MPI_DOUBLE,from,iter*10+2,MPI_COMM,&info);
  MPI_Recv(coords,NAA*35,MPI_DOUBLE,from,iter*10+3,MPI_COMM,&info);
  MPI_Recv(etcs,NAA,MPI_INT,from,iter*10+4,MPI_COMM,&info);
  MPI_Recv(nums,NAA,MPI_INT,from,iter*10+5,MPI_COMM,&info);
  MPI_Recv(ids,NAA,MPI_CHAR,from,iter*10+6,MPI_COMM,&info);
  MPI_Recv(chainids,NAA,MPI_INT,from,iter*10+7,MPI_COMM,&info);
  for(j = 0; j < NAA; j++){
    nsconformation->aa[j].etc = etcs[j];
    nsconformation->aa[j].num = nums[j];
    nsconformation->aa[j].id = ids[j];
    nsconformation->aa[j].chainid = chainids[j];
    for(i = 0; i < 3; i++){
      nsconformation->aa[j].h[i] = coords[counter++];		
	  nsconformation->aa[j].n[i] = coords[counter++];		
	  nsconformation->aa[j].ca[i] = coords[counter++];		
	  nsconformation->aa[j].c[i] = coords[counter++];		
	  nsconformation->aa[j].o[i] = coords[counter++];
	  nsconformation->aa[j].cb[i] = coords[counter++];
	  nsconformation->aa[j].g[i] = coords[counter++];
	  nsconformation->aa[j].g2[i] = coords[counter++];
	  nsconformation->xaa[j][i][0] = coords[counter++]; 
	  nsconformation->xaa[j][i][1] = coords[counter++]; 
	  nsconformation->xaa[j][i][2] = coords[counter++]; 
    }		
    nsconformation->aa[j].chi1 = coords[counter++];
    nsconformation->aa[j].chi2 = coords[counter++];
  }
  MPI_Recv(&Nchains,1,MPI_INT,from,iter*10+8,MPI_COMM,&info);
  double xaa_prev[(Nchains+1)*9];
  MPI_Recv(xaa_prev,(Nchains+1)*9,MPI_DOUBLE,from,iter*10+9,MPI_COMM,&info);
  counter = 0;
  for(j = 0; j <= Nchains; j++){
    for(i = 0; i < 3; i++){
	  nsconformation->xaa_prev[j][i][0] = xaa_prev[counter++]; 
	  nsconformation->xaa_prev[j][i][1] = xaa_prev[counter++]; 
	  nsconformation->xaa_prev[j][i][2] = xaa_prev[counter++]; 
    }
  }
  
  

}

#endif

void copyhash(ChainHash * to, ChainHash *from){
  to->index = from->index;
  to->processor = from->processor;
  to->ll = from->ll;	
}

/* A bit of ordering of the chains wrt. the logL in ascending order.
   By the end, the index with the smallest logL will be moved to the front. */
void constructhashheap(ChainHash * chainhash, int N){
  int i, booly, j;
  for(i = N/2; i > 0; i--){
    booly = 0;
    int k = i;
    ChainHash v; 
    copyhash(&v,&(chainhash[k]));
    while((booly == 0) && 2 * k <= N){
	  j = 2 * k;
	  if(j < N) //2 children
	  	if(chainhash[j].ll > chainhash[j+1].ll) j += 1;
	  if(v.ll < chainhash[j].ll) booly = 1;
	  else{
		  copyhash(&(chainhash[k]),&(chainhash[j])); 
		  k = j;
	  } 
	}	  
	copyhash(&(chainhash[k]),&v);
  }		
}
