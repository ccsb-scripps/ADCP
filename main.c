/*
** Metropolis Monte Carlo sampling procedure for simplified polypeptides.
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

#ifdef PARALLEL
#include<mpi.h>
#include"random16.h"
#endif

#include"canonicalAA.h"
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
#include"nested.h"
#include"checkpoint_io.h"
#include"flex.h"

#define VER "PepTide 2.0, Copyright (c) 2004 - 2010 Alexei Podtelezhnikov\n\
Nested Sampling by N.Burkoff\n"
#define USE "Usage: %s [options] [[-f] infile | SEQuenCE] [-o outfile]\n\
Options:\n\
 -f infile            input PDB file with initial conformation\n\
 or SEQuenCE          peptide sequence in ALPHA and beta states\n\
 -o outfile           redirected output file\n\
 -a ACCEPTANCE        crankshaft rotation acceptance rate\n\
 -A AMPLITUDE,FIX_AMP crankshaft rotation amplitude and whether the amplitude should be kept fixed (default is no fixing) \n\
 -b BETA1-BETA2:INT   thermodynamic beta schedule\n\
 -p STRING            undocumented custom parameters\n\
 -d vdw_model         Which vdW model to set the defaults for (default: LJ, could also be hard_cutoff and LJ_hard_cutoff)\n\
 -r PACExSTRETCH      test interval x total number\n\
 -s SEED              random seed\n\
 -t MASK,OPTIONS      hexadecimal mask of active tests\n\
 -c TEMP			  temperature (Celcius) to run serial MC simulation\n\
 \n\
 -n                   nested sampling procedure\n\
 -c	TEMP (NS)	  	  minimum temperature (Celcius)\n\
 -m NUM               number of MC moves between NS sample points\n\
 -r OUTPUTxMAX (NS)   output every OUTPUT sample point x max number of NS iterations\n\
 -C NUM,FILENAME      checkpointing: number NS iterations until checkpoint,checkpointfilename\n\
 -R N                 restart checkpointing - need -C, NUM=checkpointfile number\n"


#define M_LOG2E        1.4426950408889634074   


#ifdef PARALLEL
/* parallel job parameters */
int size = 1, rank = 0;

/* set thermodynamic beta for parallel tempering replicas */
void thermoset(simulation_params *sim_params)
{
	int i;
	double beta = 1.0, factor;

	if (size > 1 && sim_params->beta2 > 0.0 && sim_params->beta1 > 0.0)
		factor = pow(sim_params->beta2 / sim_params->beta1, 1.0 / (size - 1));
	else
		factor = 1.0 - 0.5 * M_LOG2E * log(size + 1) / size;

	if (sim_params->beta1 > 0.0)
		beta = sim_params->beta1;

	for (i = 0; i < rank; i++)
		beta *= factor;

	sim_params->thermobeta = beta;
	sim_params->bstp = factor;

    

	fprintf(stderr, "Rank = %d  ThermoBeta = %g\n", rank, sim_params->thermobeta);
}

void thermoswap(Chain* chain, simulation_params *sim_params)
{
	int i, j;
	static int swap = 0;
	double loss, p;
	double send[2], recv[2];
	MPI_Status status;

	swap++;

	/* who's up for a swap and what's the probability? */
	j = rand16();
	i = (j >> 8) % size;
	j = j % size;
	p = rand16();

	if ((rank != i && rank != j) || i == j)
		return;

	send[0] = sim_params->thermobeta;
	send[1] = totenergy(chain);

	if (rank == i) {
		MPI_Recv(recv, 2, MPI_DOUBLE, j, j * size + i, MPI_COMM_WORLD,
			 &status);
		MPI_Send(send, 2, MPI_DOUBLE, j, i * size + j, MPI_COMM_WORLD);
	}

	if (rank == j) {
		MPI_Send(send, 2, MPI_DOUBLE, i, j * size + i, MPI_COMM_WORLD);
		MPI_Recv(recv, 2, MPI_DOUBLE, i, i * size + j, MPI_COMM_WORLD,
			 &status);
	}

	loss = (recv[0] - send[0]) * (recv[1] - send[1]);
	if (loss < 0.0 && exp(loss) * RAND16_MAX < p)
		return;

	sim_params->thermobeta = recv[0];	/* swap */

	fprintf(stderr, "swap %d : rank %d : %g x %g <=> %g x %g\n",
		swap, rank, send[0], send[1], recv[0], recv[1]);
}

#else
void thermoset(simulation_params *sim_params)
{
	if (sim_params->beta2 > 0.0 && sim_params->beta1 > 0.0 && sim_params->pace > 0 && sim_params->stretch > 0 && sim_params->intrvl > 0)
		sim_params->bstp = pow(sim_params->beta2 / sim_params->beta1, sim_params->intrvl / ((double)sim_params->stretch));
	else
		sim_params->bstp = 1.0;

	if (sim_params->beta1 > 0.0 && sim_params->lowtemp == 0)
		sim_params->thermobeta = sim_params->beta1;
		
}

void thermoswap(Chain *chain,simulation_params *sim_params)
{
	static int swap = 0;
	//fprintf(stderr, "annealing %d : %g \n ", swap, sim_params->thermobeta);
	if (sim_params->bstp == 1.0)
		return;

	swap++;
	sim_params->thermobeta *= sim_params->bstp;

	fprintf(stderr, "annealing %d : %g\n", swap, sim_params->thermobeta);
}
#endif

static double calculateRMSD(Chain *chain, Chain *chain2)
{
	//return 3.5;
	double RMSD = 0;
	double dist = 0;
	for (int i = 0; i < chain->NAA; i++) {
		dist = (chain->aa[i].ca[0] - chain2->aa[i].ca[0])*(chain->aa[i].ca[0] - chain2->aa[i].ca[0]);
		dist += (chain->aa[i].ca[1] - chain2->aa[i].ca[1])*(chain->aa[i].ca[1] - chain2->aa[i].ca[1]);
		dist += (chain->aa[i].ca[2] - chain2->aa[i].ca[2])*(chain->aa[i].ca[2] - chain2->aa[i].ca[2]);
		RMSD += sqrt(dist);
	}
	return RMSD / (chain->NAA-1);	
}

void simulate(Chain * chain, Chaint *chaint, Biasmap* biasmap, simulation_params *sim_params)
{
	unsigned int i, j, k = sim_params->intrvl;
	double temp;
	Chain *chain2 = (Chain *)malloc(sizeof(Chain));
	chain2->aa = NULL; chain2->xaa = NULL; chain2->erg = NULL; chain2->xaa_prev = NULL;
	allocmem_chain(chain2,chain->NAA,chain->Nchains);
	if (sim_params->protein_model.external_potential_type == 5 && sim_params->protein_model.external_r0[0] == (int)sim_params->protein_model.external_r0[0]) {
		targetBest = 99999.;
		double targetBestTemp = targetBest;
        //double targetBestPrev = targetBest;
		double currTargetEnergy = 99999.;
		double lastTargetEnergy = 9999.;
		int lastIndex = 0;
		int bestIndex = 0;
		int lastBestIndex = 0;
		int pdbCount = 0;
		int stopSignal = 0;
		int swapInd = 0;
		int swapInd2 = 0;
		int swapLength = 10;
		int inCache = 0;
		int ind = 0;
		FILE *swapFile = NULL;
		char swapname[12];
		sprintf(swapname, "swap%d.pdb", swapLength);
		double swapEnergy[swapLength + 1];
		Chain *lastChain = (Chain *)malloc(sizeof(Chain));
		lastChain->aa = NULL; lastChain->xaa = NULL; lastChain->erg = NULL; lastChain->xaa_prev = NULL;

		allocmem_chain(lastChain, chain->NAA, chain->Nchains);
		Chain* swapChains[swapLength];

		for (int i = 0; i < swapLength+1; i++) {
			swapChains[i] = (Chain *)malloc(sizeof(Chain));
			swapChains[i]->aa = NULL; swapChains[i]->xaa = NULL; swapChains[i]->erg = NULL; swapChains[i]->xaa_prev = NULL;
			allocmem_chain(swapChains[i], chain->NAA, chain->Nchains);
			copybetween(swapChains[i], chain);
			swapEnergy[i] = 9999.;
		}
		double external_k = sim_params->protein_model.external_k[0];
		//fprintf(stderr, "enefffrgies C %g %g %g\n", sim_params->protein_model.external_k[0], externalBest, extenergy(chain));
		fprintf(stderr, "begin run with pace %d stretch %d \n", sim_params->pace, sim_params->stretch);
		for (i = 1; i < sim_params->stretch; i++) {
			if (stopSignal) break;
			targetBestTemp = targetBest;
			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 100 == 1 && i < 1000) || (i % 1000 == 1)) {
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					copybetween(chain2, chain);
					move(chain2, chaint, biasmap, 0.0, &temp, -1, sim_params);
					for (j = 1; (j < sim_params->pace || j < 1024); j++) {
						move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					}
				}
			}
			targetBest = targetBestTemp;
			//fprintf(stderr, "energies C %g %g %g\n", sim_params->protein_model.external_k[0], externalBest, extenergy(chain));
			for (j = 1; (j < sim_params->pace || j < 1024); j++) {
				//targetBestPrev = targetBest;
				move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				inCache = 0;
				currTargetEnergy = extenergy(chain);
				//currTargetEnergy = totenergy(chain);
				//currTargetEnergy = targetenergy(chain);

				if (currTargetEnergy - targetBest < -0.1) {
					fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					lastIndex = i*sim_params->pace + j;

					pdbCount++;
					copybetween(lastChain, chain);
					lastTargetEnergy = currTargetEnergy;
					targetBest = currTargetEnergy;
					if (sim_params->protein_model.external_k[0] != external_k) {
						fprintf(stderr, "best energy reset temp%g %g %g \n", swapEnergy[0], swapEnergy[swapLength], currTargetEnergy);
						sim_params->protein_model.external_k[0] = external_k;
					}
					bestIndex = i*sim_params->pace + j;
					lastBestIndex = i*sim_params->pace + j;
					if (rand() % 100 < 50) {
						for (ind = 0; ind < swapLength; ++ind) {
							if (swapEnergy[ind] == 9999.) swapInd = ind;
							if (calculateRMSD(swapChains[ind], chain) < 2.0) {
								inCache = 1;
								//fprintf(stderr, "inCache");
								if (currTargetEnergy < swapEnergy[ind] - 2 ) {
                                    //fprintf(stderr, "inCache %d \n",ind);
									fprintf(stderr, "swap between best curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[ind], targetBest);
									copybetween(swapChains[ind], chain);
									swapEnergy[ind] = currTargetEnergy;
									sprintf(swapname, "swap%d.pdb", ind);
									swapFile = fopen(swapname, "w+");
									pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
								}
								break;
							}
						}
						//fprintf(stderr, "swapped");
						if (!inCache) {
							if (swapEnergy[swapInd] != 9999.) {
								swapInd = rand() % swapLength;
								swapInd2 = rand() % swapLength;
								swapInd = swapEnergy[swapInd] > swapEnergy[swapInd2] ? swapInd : swapInd2;
							}
							fprintf(stderr, "swap in best curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
							sprintf(swapname,"swap%d.pdb", swapInd);
							
							swapFile = fopen(swapname, "w+");
							pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
							//fclose(swapname);
							//copybetween(chain, swapChains[swapInd]);
							copybetween(swapChains[swapInd], chain);
							swapEnergy[swapInd] = currTargetEnergy;
						}
					}
					copybetween(swapChains[swapLength], chain);
					swapEnergy[swapLength] = currTargetEnergy;
					sprintf(swapname, "swap%d.pdb", swapLength);
					swapFile = fopen(swapname, "w+");
					pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
					//fclose(swapname);
					//if (rand() * 8 < RAND_MAX) copybetween(bestChain, chain);
					//copybetween(bestChain, chain);
					//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), "Gary_Hack.pdb", totenergy(chain));
				}
				//(targetBest > 20 || currTargetEnergy - targetBest < 3.0) &&
				else if ((i*sim_params->pace + j - lastBestIndex) > 10000000) {
					fprintf(stderr, "No improvement after %d runs last best %d, stops here.\n", i*sim_params->pace + j, bestIndex);
					stopSignal = 1;
					break;
				}
				else if ((i*sim_params->pace + j - bestIndex) > 1000000 && sim_params->protein_model.external_k[0] == external_k) {
					sim_params->protein_model.external_k[0] = external_k * 0.5;
					fprintf(stderr, "No improvement after 1000000, Heat up system\n");
					swapInd = rand() % (swapLength + 1);
					while (swapEnergy[swapInd] > currTargetEnergy) swapInd = rand() % (swapLength + 1);
                    //fprintf(stderr, "swap in best %d \n", swapInd);
					//swapInd = swapLength;
					copybetween(chain, swapChains[swapInd]);
					bestIndex = i*sim_params->pace + j;
					
					fprintf(stderr, "swap out no improv curr %g swap %g best %g \n", currTargetEnergy, swapEnergy[swapInd], targetBest);
				}
				else if (currTargetEnergy != lastTargetEnergy && (i*sim_params->pace + j - lastIndex) > sim_params->pace) {
					//fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					//tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					//fprintf(stderr, "currEnergy %g best %g \n", currTargetEnergy, targetBest);
					lastIndex = i*sim_params->pace + j;
                    if (sim_params->protein_model.external_k[0] < external_k){
						sim_params->protein_model.external_k[0] = 1.02*sim_params->protein_model.external_k[0];
                    } else sim_params->protein_model.external_k[0] = external_k;
					pdbCount++;
					copybetween(lastChain, chain);
					lastTargetEnergy = currTargetEnergy;
					//if (rand() * 32 < RAND_MAX) copybetween(bestChain, chain);

					if (currTargetEnergy - targetBest <= 5.0) {
						for (ind = 0; ind < swapLength; ++ind) {
							if (swapEnergy[ind] == 9999.) swapInd = ind;
							if (calculateRMSD(swapChains[ind], chain) < 2.0) {
								inCache = 1;

								if (currTargetEnergy < swapEnergy[ind] - 2 ) {
									copybetween(swapChains[ind], chain);
									fprintf(stderr, "swap between good curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[ind], targetBest);
									swapEnergy[ind] = currTargetEnergy;
									fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
									tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
									sprintf(swapname, "swap%d.pdb", ind);
									swapFile = fopen(swapname, "w+");
									pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
									//bestIndex = i*sim_params->pace + j;
									//sim_params->protein_model.external_k[0] = external_k;
								}
								//fprintf(stderr, "inCache2");
								break;
							}
						}
						//fprintf(stderr, "swapped2");
						if (!inCache) {
							if (swapEnergy[swapInd] != 9999.) {
								swapInd = rand() % swapLength;
								swapInd2 = rand() % swapLength;
								swapInd = swapEnergy[swapInd] > swapEnergy[swapInd2] ? swapInd : swapInd2;
							}
							//copybetween(chain, swapChains[swapInd]);
							fprintf(stderr, "swap in good curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
							sprintf(swapname, "swap%d.pdb", swapInd);
							swapFile = fopen(swapname, "w+");
							pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
							//fclose(swapname);
							copybetween(swapChains[swapInd], chain);
							swapEnergy[swapInd] = currTargetEnergy;
							fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
							tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
							
							//sim_params->protein_model.external_k[0] = external_k;
							//bestIndex = i*sim_params->pace + j;
						}
						else if (rand() % 100 < 1) {
							if (1 || external_k == sim_params->protein_model.external_k[0]) {

								swapInd = rand() % (swapLength + 1);
								while (swapEnergy[swapInd] > currTargetEnergy) swapInd = rand() % (swapLength + 1);
								//swapInd = swapLength;
								fprintf(stderr, "swap out good curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
								copybetween(chain, swapChains[swapInd]);
								//sim_params->protein_model.external_k[0] = external_k;
							}
						}
					}

					if (currTargetEnergy -targetBest > 10.0 && rand() % 100 < 10) {
						if (1 || external_k == sim_params->protein_model.external_k[0]) {
							swapInd = rand() % (swapLength + 1);
							while (swapEnergy[swapInd] > currTargetEnergy) swapInd = rand() % (swapLength + 1);
							//swapInd = swapLength;
							fprintf(stderr, "swap out bad curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
							copybetween(chain, swapChains[swapInd]);
							//sim_params->protein_model.external_k[0] = external_k;
						}

					}

				}

			}
			/* annealing or tempering swap */
			if (--k == 0) {
				thermoswap(chain, sim_params);
				k = sim_params->intrvl;
				//fprintf(stderr, "swap in best %d \n", swapInd);
			}

		}
		//clean up
		freemem_chain(lastChain); free(lastChain);
		for (int i = 0; i < swapLength+1; i++) {
			freemem_chain(swapChains[i]); free(swapChains[i]);
		}

	}
	//temperature swaps
	else if (sim_params->protein_model.external_potential_type == 5 && sim_params->protein_model.external_r0[0] - (int)sim_params->protein_model.external_r0[0] < 0.2) {
		targetBest = 99999.;
		double targetBestTemp = targetBest;
		//double targetBestPrev = targetBest;
		double currTargetEnergy = 99999.;
		double lastTargetEnergy = 9999.;
		int lastIndex = 0;
		int bestIndex = 0;
		int pdbCount = 0;
		int stopSignal = 0;
		int swapInd = 0;
		int swapLength = 6;
		int currTempInd = 3;
		Chain *lastChain = (Chain *)malloc(sizeof(Chain));
		lastChain->aa = NULL; lastChain->xaa = NULL; lastChain->erg = NULL; lastChain->xaa_prev = NULL;

		allocmem_chain(lastChain, chain->NAA, chain->Nchains);
		double swapTemps[swapLength];

		for (int i = 0; i < swapLength; ++i) {
			swapTemps[i] = 0.7 + 0.1*i;
		}
		double external_k = sim_params->protein_model.external_k[0];
		fprintf(stderr, "external energy with swapping temp \n");
		fprintf(stderr, "begin run with pace %d stretch %d \n", sim_params->pace, sim_params->stretch);
		for (i = 1; i < sim_params->stretch; i++) {
			if (stopSignal) break;
			targetBestTemp = targetBest;
			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 100 == 1 && i < 1000) || (i % 1000 == 1)) {
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					copybetween(chain2, chain);
					move(chain2, chaint, biasmap, 0.0, &temp, -1, sim_params);
					for (j = 1; (j < sim_params->pace || j < 1024); j++) {
						move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					}
				}
			}
			targetBest = targetBestTemp;
			//fprintf(stderr, "energies C %g %g %g\n", sim_params->protein_model.external_k[0], externalBest, extenergy(chain));
			for (j = 1; (j < sim_params->pace || j < 1024); j++) {
				//targetBestPrev = targetBest;
				move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				//currTargetEnergy = extenergy(chain);
				currTargetEnergy = totenergy(chain);
				//currTargetEnergy = targetenergy(chain);

				if (currTargetEnergy - targetBest < 0.0) {
					fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					lastIndex = i*sim_params->pace + j;
					sim_params->protein_model.external_k[0] = external_k;
					pdbCount++;
					copybetween(lastChain, chain);
					lastTargetEnergy = currTargetEnergy;
					targetBest = currTargetEnergy;
					bestIndex = i*sim_params->pace + j;
				}
				//(targetBest > 20 || currTargetEnergy - targetBest < 3.0) &&
				else if ((i*sim_params->pace + j - bestIndex) > 4000000) {
					fprintf(stderr, "No improvement after %d runs last best %d, stops here.\n", i*sim_params->pace + j, bestIndex);
					stopSignal = 1;
					break;
				}
				//else if ((i*sim_params->pace + j - bestIndex) > 100000) {
				//	sim_params->protein_model.external_k[0] = external_k * 1.5;
				//}
				else if (currTargetEnergy != lastTargetEnergy && (i*sim_params->pace + j - lastIndex)>sim_params->pace) {
					//fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					//tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					lastIndex = i*sim_params->pace + j;

					pdbCount++;
					copybetween(lastChain, chain);
					lastTargetEnergy = currTargetEnergy;

					if (currTargetEnergy - targetBest < 5.0) {
						fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
						tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
						bestIndex = i*sim_params->pace + j;
						if (currTempInd == 0) {
							sim_params->protein_model.external_k[0] = external_k * swapTemps[1];
							currTempInd = 1;
						}
						else if (currTempInd == swapLength - 1) {
							sim_params->protein_model.external_k[0] = external_k * swapTemps[swapLength - 2];
							currTempInd = swapLength - 2;
						}
						else {
							if (rand() * 2 < RAND_MAX) {
								sim_params->protein_model.external_k[0] = swapTemps[currTempInd + 1];
								currTempInd = currTempInd +1;
							}
							else {
								sim_params->protein_model.external_k[0] = swapTemps[currTempInd - 1];
								currTempInd = currTempInd - 1;
							}

						}
					}
				}
			}
			/* annealing or tempering swap */
			if (--k == 0) {
				thermoswap(chain, sim_params);
				k = sim_params->intrvl;


			}
		}
		freemem_chain(lastChain); free(lastChain);

	}
	else if (sim_params->protein_model.external_potential_type == 5 && sim_params->protein_model.external_r0[0] - (int)sim_params->protein_model.external_r0[0] < 0.4) {
		targetBest = 99999.;
		double targetBestTemp = targetBest;
		//double targetBestPrev = targetBest;
		double currTargetEnergy = 99999.;
		double lastTargetEnergy = 9999.;
		int lastIndex = 0;
		int bestIndex = 0;
		int pdbCount = 0;
		int stopSignal = 0;

		Chain *lastChain = (Chain *)malloc(sizeof(Chain));
		lastChain->aa = NULL; lastChain->xaa = NULL; lastChain->erg = NULL; lastChain->xaa_prev = NULL;

		allocmem_chain(lastChain, chain->NAA, chain->Nchains);

		double external_k = sim_params->protein_model.external_k[0];
		fprintf(stderr, "external energy with swapping temp \n");
		fprintf(stderr, "begin run with pace %d stretch %d \n", sim_params->pace, sim_params->stretch);
		for (i = 1; i < sim_params->stretch; i++) {
			if (stopSignal) break;
			targetBestTemp = targetBest;
			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 100 == 1 && i < 1000) || (i % 1000 == 1)) {
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					copybetween(chain2, chain);
					move(chain2, chaint, biasmap, 0.0, &temp, -1, sim_params);
					for (j = 1; (j < sim_params->pace || j < 1024); j++) {
						move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					}
				}
			}
			targetBest = targetBestTemp;
			//fprintf(stderr, "energies C %g %g %g\n", sim_params->protein_model.external_k[0], externalBest, extenergy(chain));
			for (j = 1; (j < sim_params->pace || j < 1024); j++) {
				//targetBestPrev = targetBest;
				move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				//currTargetEnergy = extenergy(chain);
				currTargetEnergy = totenergy(chain);
				//currTargetEnergy = targetenergy(chain);

				if (currTargetEnergy - targetBest < 0.0) {
					fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					lastIndex = i*sim_params->pace + j;
					sim_params->protein_model.external_k[0] = external_k;
					pdbCount++;
					copybetween(lastChain, chain);
					lastTargetEnergy = currTargetEnergy;
					targetBest = currTargetEnergy;
					bestIndex = i*sim_params->pace + j;
					lastIndex = i*sim_params->pace + j;
				}
				//(targetBest > 20 || currTargetEnergy - targetBest < 3.0) &&
				else if ((i*sim_params->pace + j - bestIndex) > 4000000) {
					fprintf(stderr, "No improvement after %d runs last best %d, stops here.\n", i*sim_params->pace + j, bestIndex);
					stopSignal = 1;
					break;
				}
				else if ((i*sim_params->pace + j - bestIndex) > 100000) {
					sim_params->protein_model.external_k[0] = external_k * 1.5;
				}
				else if (currTargetEnergy != lastTargetEnergy && (i*sim_params->pace + j - lastIndex) > sim_params->pace) {
					//fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					//tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					lastIndex = i*sim_params->pace + j;

					pdbCount++;
					
					lastTargetEnergy = currTargetEnergy;

					if (currTargetEnergy - targetBest < 10.0) {
						fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
						tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
						lastIndex = i*sim_params->pace + j;
						copybetween(lastChain, chain);
					}
				} else if (currTargetEnergy != lastTargetEnergy && (i*sim_params->pace + j - lastIndex) > 200000) {
					copybetween(chain, lastChain);
				}
			}
			/* annealing or tempering swap */
			if (--k == 0) {
				thermoswap(chain, sim_params);
				k = sim_params->intrvl;


			}
		}
		freemem_chain(lastChain); free(lastChain);

	} else if (sim_params->protein_model.external_potential_type == 5) {
		//double targetBestPrev = targetBest;
		double targetBestTemp = targetBest;
		targetBest = 9999.;
		double currTargetEnergy = 99999.;
		double lastTargetEnergy = 9999.;
		for (i = 1; i < sim_params->stretch; i++) {
			targetBestTemp = targetBest;
			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 100 == 1 && i < 1000) || (i % 1000 == 1)) {
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					copybetween(chain2, chain);
					move(chain2, chaint, biasmap, 0.0, &temp, -1, sim_params);
					for (j = 1; (j < sim_params->pace || j < 1024); j++) {
						move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					}
				}
			}
			targetBest = targetBestTemp;
			for (j = 1; (j < sim_params->pace || j < 1024); j++) {
				//targetBestPrev = targetBest;
				move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				currTargetEnergy = extenergy(chain);
				//currTargetEnergy = targetenergy(chain);
				if (currTargetEnergy - targetBest < 0.0) {
					//fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					//tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					//lastTargetEnergy = currTargetEnergy;
					targetBest = currTargetEnergy;
					//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), "Gary_Hack.pdb", totenergy(chain));
				}
				

			}
#ifdef PARALLEL
			/*static FILE* fptr;
			if(fptr == NULL){
			char filename[9];
			sprintf(filename, "%d.temp", rank);
			fptr = fopen(filename,"w");
			//fprintf(fptr,"#T: %f\n", 1000.0/(1.9858775*thermobeta)-273.15);
			}
			fprintf(fptr,"%f %f\n",thermobeta,totenergy()); */
			if (sim_params->thermobeta != sim_params->beta1)
				continue;
#endif
			if (1) {
				fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
				tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
				//lastTargetEnergy = currTargetEnergy;
			}
			/* annealing or tempering swap */
			if (--k == 0) {
				thermoswap(chain, sim_params);
				k = sim_params->intrvl;
			}
		}
	}
	else {
		for (i = 1; i < sim_params->stretch; i++) {

			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 100 == 1 && i < 1000) || (i % 1000 == 1)) {
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					copybetween(chain2, chain);
					move(chain2, chaint, biasmap, 0.0, &temp, -1, sim_params);
					for (j = 1; (j < sim_params->pace || j < 1024); j++) {
						move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					}
				}
			}
			for (j = 1; (j < sim_params->pace || j < 1024); j++) {
				move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				/* annealing or tempering swap */
				if (--k == 0) {
					thermoswap(chain, sim_params);
					k = sim_params->intrvl;

				}
			}
#ifdef PARALLEL
			/*static FILE* fptr;
			if(fptr == NULL){
			char filename[9];
			sprintf(filename, "%d.temp", rank);
			fptr = fopen(filename,"w");
			//fprintf(fptr,"#T: %f\n", 1000.0/(1.9858775*thermobeta)-273.15);
			}
			fprintf(fptr,"%f %f\n",thermobeta,totenergy()); */
			if (sim_params->thermobeta != sim_params->beta1)
				continue;
#endif

			fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
			tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
		}
	}

	freemem_chain(chain2); free(chain2);

}

char *read_options(int argc, char *argv[], simulation_params *sim_params)
{
	unsigned int pace = 0, stretch = 16, tmask = 0x0, seed = 0;
	double thermobeta;
	int lowtemp;
	double beta1 = 1.0, beta2 = 0.0; //, bstp = 1.0;
	unsigned int intrvl = 16384;
	int NS = 0; //1 = yes
	int iter_max = 1000;
	int num_NS_per_checkpoint = 0;
	char checkpoint_filename[256];
	int checkpoint_counter = 0;
	int checkpoint = 0;
	int restart_from_checkpoint = 0;

	int i, opt;
	char *retval = NULL;
	double acceptance_rate = sim_params->acceptance_rate;
	double amplitude = sim_params->amplitude;
	int keep_amplitude_fixed = sim_params->keep_amplitude_fixed;
	char error_string[DEFAULT_LONG_STRING_LENGTH]="";


	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			//if (freopen(argv[i], "r", stdin) == NULL)
			//try to open file
			if ((sim_params->infile = fopen(argv[i],"r")) == NULL)
				//failed: it must be a sequence
				retval = argv[i];
			else
				if (sim_params->infile_name) free(sim_params->infile_name);
				copy_string(&(sim_params->infile_name),argv[i]);
			continue;
		}

		opt = argv[i][1];
		if (++i >= argc && opt != 'n')
			opt = 0;

		switch (opt) {
		case 'a':
			sscanf(argv[i], "%lf", &acceptance_rate);
			if(acceptance_rate > 1.0 || acceptance_rate <= 0.0){
			  fprintf(stderr,"Acceptance rate must be between 0.0 and 1.0, setting it to 0.5\n");
			  acceptance_rate = 0.5; 	
			}
			sim_params->acceptance_rate = acceptance_rate;			
			
			break;
		case 'A': //amplitude
			sscanf(argv[i], "%lf,%d", &amplitude,&keep_amplitude_fixed);
			if(amplitude > 0.0 || amplitude <= -M_PI){
			  fprintf(stderr,"Amplitude must be between -PI and 0.0, setting it to -0.25\n");
			  amplitude = -0.25;
			}
			sim_params->amplitude = amplitude;
			sim_params->keep_amplitude_fixed = keep_amplitude_fixed;
			break;
		case 'b':
			sscanf(argv[i], "%lf-%lf:%u", &beta1, &beta2, &intrvl);
			sim_params->beta1 = beta1;
			sim_params->beta2 = beta2;
			sim_params->intrvl = intrvl;
			break;
		case 'c':
			sscanf(argv[i],"%lf",&thermobeta);
		    thermobeta = 298.0/(thermobeta+273.0);
		    lowtemp = 1;
			sim_params->thermobeta = thermobeta;
			sim_params->lowtemp = lowtemp;
			break;
		case 'C':
			sscanf(argv[i], "%d,%256s",&num_NS_per_checkpoint,checkpoint_filename);
			checkpoint = 1;
			sim_params->num_NS_per_checkpoint = num_NS_per_checkpoint;
			sim_params->checkpoint_filename = realloc(sim_params->checkpoint_filename,DEFAULT_SHORT_STRING_LENGTH);
			strcpy(sim_params->checkpoint_filename,checkpoint_filename);
			sim_params->checkpoint = checkpoint;
			break;
		case 'f':
			//if (freopen(argv[i], "r", stdin) == NULL)
			//try to open file
			if ((sim_params->infile = fopen(argv[i],"r")) == NULL)
				//failed: it must be a sequence
				retval = argv[i];
			else
				if (sim_params->infile_name) free(sim_params->infile_name);
				copy_string(&(sim_params->infile_name),argv[i]);
			break;
		case 'm':
		    sscanf(argv[i], "%u", &iter_max);
			sim_params->iter_max = iter_max;
			break;
		case 'M':
		    sscanf(argv[i],"%u",&(sim_params->number_initial_MC));
		    break;
		case 'n':
		    NS = 1;
			sim_params->NS = NS;
		    i--;
		    break;
		case 'o':
			//freopen(argv[i], "w", stdout);
			if ((sim_params->outfile = fopen(argv[i],"w")) == NULL)
				stop("Could not open output file.\n");
			else
				if (sim_params->outfile_name) free(sim_params->outfile_name);
				copy_string(&(sim_params->outfile_name),argv[i]);
			break;
		case 'p':
			if (sim_params->prm) free(sim_params->prm);
			copy_string(&sim_params->prm,argv[i]);
			break;
		case 'd':
			if (strcmp(argv[i],"hard_cutoff")==0) {
				sim_params->protein_model.vdw_potential=HARD_CUTOFF_VDW_POTENTIAL;
				set_hard_cutoff_default_params(&(sim_params->protein_model));
			} else if (strcmp(argv[i],"lj")==0) {
				sim_params->protein_model.vdw_potential=LJ_VDW_POTENTIAL;
				set_lj_default_params(&(sim_params->protein_model));
			} else if (strcmp(argv[i],"lj_hard_cutoff")==0) {
				sim_params->protein_model.vdw_potential=LJ_VDW_POTENTIAL;
				set_lj_default_params(&(sim_params->protein_model));
				sim_params->protein_model.vdw_lj_neighbour_hard=1;
				sim_params->protein_model.vdw_lj_hbonded_hard=1;
			} else {
				sprintf(error_string,"Unknown value for vdW potential model (%s).  It must be one of hard_cutoff and lj.",argv[i]);
				stop(error_string);
			}
			break;
		case 'r':
			sscanf(argv[i], "%ux%u", &pace, &stretch);
			sim_params->pace = pace;
			sim_params->stretch = stretch;
			break;
		case 'R':
			sscanf(argv[i], "%d",&checkpoint_counter);
			if(checkpoint_counter != -1){
			  restart_from_checkpoint = 1;
			  sim_params->checkpoint_counter = checkpoint_counter;
			  sim_params->restart_from_checkpoint = restart_from_checkpoint;
			}
			break;
		case 's':
			sscanf(argv[i], "%u", &seed);
			sim_params->seed = seed;
			break;
		case 't':
			sscanf(argv[i], "%x", &tmask);
			sim_params->tmask = tmask;
			break;
		default:
			fprintf(stderr, VER USE PARAM_USE, argv[0]);
			helps();
			exit(EXIT_FAILURE);
		}
	}

	if (sim_params->seq) free(sim_params->seq);
	copy_string(&(sim_params->seq),retval);

	return sim_params->seq;
}

void graceful_exit(int sig)
{
	exit(EXIT_FAILURE);	/* flushes streams too */
}

void single_point_test(Chain *chain, Chaint *chaint, Biasmap *biasmap, simulation_params *sim_params,unsigned int i)
{

	//fprintf(stderr,"SINGLE POINT\n");
	//fprintf(stderr,"XAA\n");
	//for (int i = 1; i < chain->NAA; i++) {
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");

	/* correct peptide if needed */
	if(sim_params->protein_model.fixit) fixpeptide(chain->aa, chain->NAA, &(sim_params->protein_model)); // PEPTIDE MODIFICATION!!
	chkpeptide(chain->aa, chain->NAA, &(sim_params->protein_model));
	update_sim_params_from_chain(chain,sim_params); // updating NAA and seq

	biasmap_initialise(chain,biasmap,&(sim_params->protein_model));
	energy_matrix_calculate(chain,biasmap,&(sim_params->protein_model));

	/* tests before projecting the peptide onto the CRANKITE model */
	fprintf(sim_params->outfile,"-+- PLAY BLOCK %5d -+-\n", i);
//	fprintf(sim_params->outfile,"-+- before init %5d -+-\n", i);
	tests(chain,biasmap,sim_params->tmask, sim_params, 0x1, NULL );
//	fprintf(sim_params->outfile,"-+- end before init %5d -+-\n", i);


	initialize(chain,chaint,sim_params); // PEPTIDE MODIFICATION!!

	/* tests before projecting the peptide onto the CRANKITE model */
	energy_matrix_calculate(chain,biasmap,&(sim_params->protein_model));
//	fprintf(sim_params->outfile,"-+- after init %5d -+-\n", i);
	tests(chain,biasmap,sim_params->tmask, sim_params, 0x10, NULL );
//	fprintf(sim_params->outfile,"-+- end after init %5d -+-\n", i);

}

void set_random_seed(simulation_params *sim_params) {

	if (sim_params->seed == 0)
		sim_params->seed = (unsigned int) time(NULL);

	/*random seed*/
#ifdef PARALLEL
	srand(sim_params->seed + 257 * rank);	/* asynchronous RNG */
	/* if MC, not NS, broadcast synchronous RNG */
	if(!sim_params->NS){
	MPI_Bcast(&(sim_params->seed), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	  srand16(sim_params->seed);		/* synchronous RNG */
	}
#else
	srand(sim_params->seed);
#endif

}

int main(int argc, char *argv[])
{
	//char *seq;
	simulation_params sim_params;

	
	signal(SIGTERM, graceful_exit);

#ifdef PARALLEL
	/* initialise MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	
	/* SET SIMULATION PARAMS */
	param_initialise(&sim_params); //set default
	set_lj_default_params(&(sim_params.protein_model)); // set the default parameters for the LJ model
	read_options(argc, argv, &sim_params);

	//param_print(sim_params,sim_params.outfile); //default + read-in

	set_random_seed(&sim_params);

	//set a timer
	time_t startTime = time(NULL);

	/* set different default simulation params for NS and MC */
	if(sim_params.NS){
	  if(sim_params.tmask == 0x0) sim_params.tmask = 0x18803;
	  if(sim_params.pace == 0) {sim_params.pace = 1; sim_params.stretch = 1000;}
	} else {
	  if(sim_params.tmask == 0x0) sim_params.tmask = 0x3;
	  sim_params.tmask |= 0x8800;
	  if(sim_params.pace==0) {sim_params.pace = 4096; sim_params.stretch = 16;} 
	}

	/* set temperature for MC */
	if(!sim_params.NS){
	  thermoset(&sim_params);
	}

	/* INITIALISE PROTEIN MODEL */

	//model_param_initialise(&(sim_params.protein_model));
	model_param_read(sim_params.prm,&(sim_params.protein_model),&(sim_params.flex_params));



	initialize_sidechain_properties(&(sim_params.protein_model));
	vdw_cutoff_distances_calculate(&sim_params, stderr, 0);
	peptide_init();



	param_print(sim_params,sim_params.outfile); //read-in
    
	if (sim_params.protein_model.external_potential_type == 5) {

		Cgridmap_initialise();
		fprintf(stderr, "C Grid map initialisation finished");
		fprintf(stderr, "best target energy %g\n", targetBest);
		CAgridmap_initialise();
		fprintf(stderr, "CA Grid map initialisation finished");
		NAgridmap_initialise();
		Hgridmap_initialise();
		Sgridmap_initialise();
		Ogridmap_initialise();
		//fprintf(stderr, "O Grid map initialisation finished");
		Ngridmap_initialise();
		egridmap_initialise();
		dgridmap_initialise();
		fprintf(stderr, "AD Grid maps initialisation finished \n");
	}
	/* HERE STARTS THE ACTUAL SIMULATION */

	if(!sim_params.NS){
	    /* allocate memory for the peptide */
	    Chain *chain = (Chain *)malloc(sizeof(Chain)); chain->NAA = 0;
            Chaint* chaint = (Chaint *)malloc(sizeof(Chaint));
      	    chain->NAA =0;
      	    chain->aa = NULL; chain->xaa = NULL; chain->erg = NULL; chain->xaa_prev = NULL;
      	    chaint->aat = NULL; chaint->xaat = NULL; chaint->ergt = NULL; chaint->xaat_prev = NULL;
	    /* allocate memory for the biasmap */
            Biasmap *biasmap = (Biasmap *)malloc(sizeof(Biasmap));
      	    biasmap->distb = NULL;

	   /* read in / generate the peptide */
	   if (sim_params.seq != NULL) {

		//if (sim_params.protein_model.external_constrained_aalist_file) {
		//	stop("constraints unimplemented!");
		//}
		/* build peptide from scratch, do not do tests */
		build_peptide_from_sequence(chain,chaint,sim_params.seq, &sim_params);
		mark_fixed_aa_from_file(chain,&sim_params);
		mark_constrained_aa_from_file(chain,&sim_params);
		chkpeptide(chain->aa, chain->NAA, &(sim_params.protein_model));
		update_sim_params_from_chain(chain,&sim_params); // updating NAA and seq
		biasmap_initialise(chain,biasmap,&(sim_params.protein_model));
		energy_matrix_calculate(chain,biasmap,&(sim_params.protein_model));
	    
	    
	   } else { /* read in peptide */

		unsigned int i;
		int readin = 0;

		if (sim_params.stretch == 0) { /* single point testing */

#ifdef PARALLEL
		   if (rank==0) {  // only do the tests and printing on the master node
		   fprintf(stderr,"WARNING! Single point tests are done on the master node only!  A parallel call with stretch=0 might waste resources.\n");
#endif
			fprintf(stderr,"INFO: 0 MC step was asked for.  Only doing tests on the PDB entries.\n");
			/* while reading, do test one by one, to save memory */

			for (i = 1; pdbin(chain,&sim_params,sim_params.infile) != EOF; i++) {
			    mark_fixed_aa_from_file(chain,&sim_params);
			    mark_constrained_aa_from_file(chain,&sim_params);
			    single_point_test(chain,chaint,biasmap,&sim_params,i);
			}
#ifdef PARALLEL
		    }
#endif

		} else { /* reading in for an MC simulation */
#ifdef PARALLEL
			fprintf(stderr,"INFO: Attempting a parallel tempering simulation.\n");
			/* parallel tempering: upto max(rank) initial configs */
			int i;
			for (i = 0; pdbin(chain,&sim_params,sim_params.infile) != EOF && i < rank; i++)  readin = 1;
			if (readin == 0)  fprintf(stderr,"WARNING! Rank %d did not read in a PDB\n",rank);
        
#else
			/* serial MC:  only last entry */
			fprintf(stderr,"INFO: Attempting a serial MC simulation.\n");
			unsigned int i;

			for (i = 1; pdbin(chain,&sim_params,sim_params.infile) != EOF; i++);
			if (i>2) fprintf(stderr,"WARNING! First %d entries in the input PDB will be ignored, only using last one for the MC simulation.\n",i-1);
			if (i==1) {
				stop("ERROR! EOF while reading in from input PDB file.");
			}
			readin=1;
#endif

		/* project the peptide onto the CRANKITE model and do the initial tests for serial MC */
			if (readin) {
				if(sim_params.protein_model.fixit) fixpeptide(chain->aa, chain->NAA, &(sim_params.protein_model));
				chkpeptide(chain->aa, chain->NAA, &(sim_params.protein_model));
				mark_fixed_aa_from_file(chain,&sim_params);
				mark_constrained_aa_from_file(chain,&sim_params);
			}
			update_sim_params_from_chain(chain,&sim_params); // updating NAA and seq
			//fprintf(stderr,"Updated sim_params->seq: %s\n",sim_params.seq); //this also has a starting A which is not part of the polypeptide
#ifndef PARALLEL
			fprintf(sim_params.outfile,"-+- PLAY BLOCK %5d -+-\n", --i);
			biasmap_initialise(chain,biasmap,&(sim_params.protein_model));
			energy_matrix_calculate(chain,biasmap,&(sim_params.protein_model));
			tests(chain,biasmap,sim_params.tmask, &sim_params, 0x1, NULL );
#endif
			initialize(chain,chaint,&sim_params); // peptide modification
			/* allocate energy matrix and read in biasmap */
			biasmap_initialise(chain,biasmap,&(sim_params.protein_model));
			energy_matrix_calculate(chain,biasmap,&(sim_params.protein_model));
#ifndef PARALLEL
			/* initial test */
			tests(chain,biasmap,sim_params.tmask, &sim_params, 0x10, NULL );
#endif
		}
	  }



	/* MC */
	simulate(chain,chaint,biasmap,&sim_params);
	finalize(chain,chaint,biasmap); //free memory allocated in initialize
     
	} else { /* Nested Sampling, parallel or serial */
	  nestedsampling(sim_params.pace,sim_params.stretch,&sim_params);
	}
	
	param_finalise(&sim_params);


#ifdef PARALLEL
    /* finalise MPI */
    MPI_Finalize();
	if (rank == 0)
#endif
	free(Cmapvalue);
	free(CAmapvalue);
	free(NAmapvalue);
	free(Smapvalue);
	free(Hmapvalue);
	free(Omapvalue);
	free(Nmapvalue);
	free(emapvalue);
	free(dmapvalue);
	fprintf(stderr, "best target energy %g\n", targetBest);


	fprintf(stderr,"The program has successfully finished in %d seconds. :)  Bye-bye!\n", time(NULL)- startTime);
	return EXIT_SUCCESS;
}
