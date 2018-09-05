/* Nested Sampling procedure for Crankite
 * 
 * Skeleton of code from J.Skilling p188 of
 * Data Analysis A Bayesian Tutorial D.S.Sivia
 *
 * Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
 * 
 */

typedef struct{
  int *instructions;
  int length;
  int current_position;
} Instructions;


/*Instruction Help
 *
 * Example 1:
 *
 * -1 5 743  //receive chain 743 from proc 5
 *  0 0 -1   //run MC moves
 *  0 0 23   //put MC'd chain in position 23 on this processor
 *
 *
 * Example 2:
 *
 * 1 6 35 //send chain 35 to proc 6
 * 1 4 98 //send chain 98 to proc 4
 * 0 0 46 //copy chain 46 to run MC
 * 1 3 78 // send MC'd chain into proc 3 index  78
 *
 *
 * Example 3:
 *  0 0 46 //copy chain 46 to run MC
 *  0 0 35 //put MC'd chain in position 35 on this processor
 *
 * Example 4:
 *  1 8 100 //receive chain from proc 8 position 100 to MC
 *  0 0 -1 // MC chain
 *  -1 4 67 // receive MC'd chain from proc 4 put in index 67
 *  -1 9 111 // receive MC'd chain from proc 9 put in index 111
 *  0 0 56 // put local MC'd chain in index 56 (on this proc)
 */

void nestedsampling(int,int, simulation_params *);

