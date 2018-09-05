/*
** Module of error handling **
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdio.h>
#include<stdlib.h>
//#include<string.h>
#ifdef PARALLEL
#include<mpi.h>
#endif
#include"error.h"

/* Aborting the program with an error message   */
/* Different failure codes could be passed here */
//void stop(char *error_string, int failure_code) {
void stop(char *error_string) {

   int failure_code = EXIT_FAILURE;

#ifdef PARALLEL

   int myrank;
   MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

   fprintf(stderr,"ERROR! (MPI rank %d) %s\n",myrank,error_string);
   MPI_Abort(MPI_COMM_WORLD, failure_code);

#else

   fprintf(stderr,"ERROR! %s\n",error_string);
   exit(failure_code);

#endif

}
