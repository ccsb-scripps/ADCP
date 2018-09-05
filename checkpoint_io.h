/*
** Checkpoint input-output routines.
**
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

typedef struct _ChainHash{
  int processor;
  int index;
  double ll;	
} ChainHash;
int read_in_from_checkpoint(simulation_params *sim_params,
			Biasmap **biasmap,
			Chain *temporary,
			Chain **cpoints,
			Chaint **chaint,
			ChainHash **chainhash,
			int rank,
			int P,void*mpi_comm,int only_output_checkpoint);
int read_in_from_pdb(simulation_params *sim_params,
			Biasmap **biasmap, Chain *temporary, Chain **cpoints, Chaint **chaint,
			ChainHash **chainhash, int rank, int P, void *mpi_comm);
void initialize_all_pdb_chains(simulation_params *sim_params,
			Biasmap **biasmap, Chain *temporary, Chain **cpoints, Chaint **chaint,
			ChainHash **chainhash, int rank, int current_stored, int P, void *mpi_comm);
void open_next_checkpoint_file(simulation_params *sim_params, char *action);
void initialize_checkpoint_file_as_file_pointer_to_print(simulation_params *sim_params, char *action);
void read_checkpoint_header(simulation_params *sim_params);
void read_checkpoint_entry(Chain *cpoints, simulation_params *sim_params);
void print_checkpoint_header(simulation_params *sim_params, FILE * outfile, int j);
void print_checkpoint_entry(Chain *cpoints, simulation_params *sim_params, FILE *outfile, int N);
void output_checkpoint_file(Chain *cpoints, int current_stored, simulation_params *sim_params,void *mpi_comm);
#ifdef PARALLEL
void mpi_send_chain(Chain* nsconformation, int from, int to, double *logLstar, int iter,MPI_Comm MPI_COMM);
void mpi_rec_chain(Chain *nsconformation, int from, int to, double *logLstar, int iter,MPI_Comm MPI_COMM);
#endif
void copyhash(ChainHash * to, ChainHash *from);
void constructhashheap(ChainHash * chainhash, int N);
int store_chain(ChainHash **chainhash, Chain *temporary, Biasmap *biasmap, Chaint *chaint, int P, int *current_stored, int counter, Chain **cpoints, simulation_params *sim_params, int rank, void*mpi_comm);
