/*
** This is a core of Metropolis Monte Carlo sampling procedure for
** simplified polypeptides. It contains the simulation procedure
** as well as the setup of initial conformation based on provided PDB
** or simplified sequence-structure input.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>		/* rand() */
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<float.h>

#include"error.h"
#include"params.h"
#include"aadict.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"
#include"vdw.h"
#include"energy.h"
#include"metropolis.h"


#define Erg(I,J)     erg[(I) * chain->NAA + (J)]
#define Ergt(I,J)   ergt[(I - start) * chain->NAA + (J)]
//#define Ergt(I,J)   ergt[(I) * chain->NAA + (J)]


/***********************************************************/
/****           MOVES AND METROPOLIS CRITERIA           ****/
/***********************************************************/

/* Modulo a number to range(1..mod)
reModNum(6,6)=6;
reModNum(1,6)=1
reModNum(7,6)=6;
reModNum(0,6)=6
*/
static int reModNum(int a, int mod) {
	return ((a-1)%mod)+1;
}

static int indMoved(int ind, int start, int end){
	if (start>=end){
		if (ind>end && ind<start)
			return 0;
		else
			return 1;
	} else {
		if (ind>=start && ind<=end)
			return 1;
		else
			return 0;
	}
}

/* Check if the proposed move (saved in chaint) is allowed
   by applying the Metropolis criteria on the energy change.
   If the move is allowed, update the coordinates and the
   energy matrix. */
static int allowed(Chain *chain, Chaint *chaint, Biasmap* biasmap, int start, int end, double logLstar, double *currE, simulation_params *sim_params)
{	
	int i, j;
	double q, loss = 0.0;
	int linked = 0;
	//get the AD energy first as it will set position for gamma atoms
	double externalloss = 0.0;
	//double ADEnergy_Chaint[end-start+1];
	double* ADEnergy_Chaint;
	if (sim_params->protein_model.external_potential_type == 5){
		
		ADEnergy_Chaint = ADenergyNoClash(start,end,chain,chaint,&(sim_params->protein_model), 0);
		for (i = start; i <= end; i++){
			//ADEnergy_Chaint[i-start] = chaint->Ergt(0, i);
			externalloss += chain->Erg(0, reModNum(i, chain->NAA-1)) - ADEnergy_Chaint[i-start];
		}	
	}
	if (end>chain->NAA-1) linked=1;
	for (i = start; i <= end; i++){
		for (j = 1; j < chain->NAA; j++) {
			if (j == reModNum(i, chain->NAA-1)){
				q = energy1(chaint->aat + j, &(sim_params->protein_model));
				if ( ( j==1 || j==chain->NAA-1)){
					if (sim_params->protein_model.external_potential_type2 != 4)
						q += 0;
					else if (j==1 && indMoved(2,start,reModNum(end,chain->NAA-1)) && linked)
						q += ramabias(chaint->aat + chain->NAA - 1, chaint->aat + 1, chaint->aat + 2);
					else if (j==1 && linked)
						q += ramabias(chaint->aat + chain->NAA - 1, chaint->aat + 1, chain->aa + 2);	
					else if (j==1)
						q += ramabias(chain->aa + chain->NAA - 1, chaint->aat + 1, chaint->aat + 2);	
					else if (j==chain->NAA-1 && indMoved(chain->NAA-1,start,reModNum(end,chain->NAA-1)) && linked)
						q += ramabias(chaint->aat + j - 1, chaint->aat + j, chaint->aat +1);
					else if (j==chain->NAA-1 && linked)
						q += ramabias(chain->aa + j - 1, chaint->aat + j, chaint->aat +1);
					else
						q += ramabias(chaint->aat + j - 1, chaint->aat + j, chain->aa +1);
				} else if (i == start)
					q += ramabias(chain->aa + reModNum(i-1, chain->NAA-1), chaint->aat + reModNum(i, chain->NAA-1), chaint->aat + reModNum(i+1, chain->NAA-1));
				else if (i == end)
					q += ramabias(chaint->aat + reModNum(i-1, chain->NAA-1), chaint->aat + reModNum(i, chain->NAA-1), chain->aa + reModNum(i+1, chain->NAA-1));
				else
					q += ramabias(chaint->aat + reModNum(i-1, chain->NAA-1), chaint->aat + reModNum(i, chain->NAA-1), chaint->aat + reModNum(i+1, chain->NAA-1));
			} 
			else if (indMoved(j,start,reModNum(end,chain->NAA-1))){
				if ((reModNum(i, chain->NAA-1) == 1 && j == chain->NAA-1)) {
					q = energy2cyclic(biasmap,chaint->aat + 1, chaint->aat + chain->NAA - 1, &(sim_params->protein_model));
					chaint->Ergt(j, reModNum(i, chain->NAA-1)) = q;
				} else if(j > reModNum(i, chain->NAA-1)) {
					//fprintf(stderr,"MC move q = %d %d, loss = %d %d haha %d %d,",i,j,start,end,indMoved(j,start,reModNum(end,chain->NAA-1)),linked);
					q = energy2(biasmap,(chaint->aat) + reModNum(i, chain->NAA-1), (chaint->aat) + j, &(sim_params->protein_model));
					if (j < start && linked) {
						//fprintf(stderr,"aaa move q = %d %d, loss = %d %d haha %d %d,\n",i,j,start,end,indMoved(j,start,reModNum(end,chain->NAA-1)),linked);
						chaint->Ergt(j + chain->NAA-1, reModNum(i, chain->NAA-1)) = q;
					} else {
						//fprintf(stderr,"MC move q = %d %d, loss = %d %d haha %d %d,\n",i,j,start,end,indMoved(j,start,reModNum(end,chain->NAA-1)),linked);
						chaint->Ergt(j, reModNum(i, chain->NAA-1)) = q;
					}
						
				} else {
					//q = energy2(biasmap,(chaint->aat) + reModNum(i, chain->NAA-1), (chaint->aat) + j, &(sim_params->protein_model));
					//chaint->Ergt(reModNum(i, chain->NAA-1), j) = chaint->Ergt(j, reModNum(i, chain->NAA-1));
					continue;
				}
			} else {
				if (j == 1 && i == chain->NAA-1)
					q = energy2cyclic(biasmap,chain->aa + 1, chaint->aat + chain->NAA - 1, &(sim_params->protein_model));
				else if (i == 1 && j == chain->NAA-1)
					q = energy2cyclic(biasmap,chaint->aat + 1, chain->aa + chain->NAA - 1, &(sim_params->protein_model));
				else
					q = energy2(biasmap,chaint->aat + reModNum(i, chain->NAA-1), (chain->aa) + j, &(sim_params->protein_model));
			}

			chaint->Ergt(i, j) = q;
			loss += (chain->Erg(reModNum(i, chain->NAA-1), j) - q);
			//fprintf(stderr,"MC move i = %d, j = %d, %g\n",i,j,loss);
		}
	}
	/*Also take into account the global_energy term */


	//q = global_energy(start,end,chain,chaint,biasmap,&(sim_params->protein_model));
	//loss += (chain->Erg(0, 0) - q);
	//fprintf(stderr,"MC move q = %g, loss = %g,",q,loss);
	//double externalloss = (chain->Erg(0, 0) - q);





	double internalloss = loss;
	double external_k = 1.0;
	if (sim_params->protein_model.external_potential_type == 5 || sim_params->protein_model.external_potential_type2 == 4)	external_k = sim_params->protein_model.external_k[0];
	//if (chain->Erg(0, 0) > 5 || currTargetEnergy - targetBest > 15) external_k = 0.5;
	//if (externalloss < -10 || loss < -10) external_k = 0.05 * external_k;
	//if (chain->Erg(0, 0) > 20 ||chain->Erg(0, 0) > 50) external_k = 0.2 * external_k;
	//loss is negative!! if loss is negative, it's worse, bad
	/* Metropolis criteria */
	//loss += q - chain->Erg(0, 0);
	//loss = loss/sqrt(chain->NAA) + externalloss;

	double cyclicBondEnergy;
	/*special cyclic*/  // needs attention !!! sign might be wrong...!!!GARY HACK
	if (sim_params->protein_model.external_potential_type2 == 4) {
		if (linked)
			cyclicBondEnergy = cyclic_energy((chaint->aat) + 1, (chaint->aat) + chain->NAA - 1, 0);
		else if (start == 1)
			cyclicBondEnergy = cyclic_energy((chaint->aat) + 1, (chain->aa) + chain->NAA - 1, 0);
		else if (end == chain->NAA-1) 
			cyclicBondEnergy = cyclic_energy((chain->aa) + 1, (chaint->aat) + chain->NAA - 1, 0);
		else
			cyclicBondEnergy = cyclic_energy((chain->aa) + 1, (chain->aa) + chain->NAA - 1, 0);
		externalloss += chain->Erg(1, 0) - cyclicBondEnergy;
	}



	//loss = loss + externalloss;


	
	

	if (loss < 0.0  && !sim_params->NS &&  exp(sim_params->thermobeta * (loss)) * RAND_MAX < rand()) {
		//fprintf(stderr," rejected\n");
		if (sim_params->protein_model.external_potential_type == 5)
			free(ADEnergy_Chaint);
		return 0;	/* disregard rejected changes */
	}
	//if (loss < 0.0 ) {
	////	//fprintf(stderr," rejected\n");
	//	return 0;	/* disregard rejected changes */
	//}

	

	//if (loss < 0.0  && !sim_params->NS &&  exp(sim_params->thermobeta * loss) * RAND_MAX < rand()) {
	//	//fprintf(stderr," rejected\n", );
	//	return 0;	/* disregard rejected changes */
	//}

	
	//if (loss < 0.0  && !sim_params->NS &&  exp(sim_params->thermobeta * externalloss *external_k) * RAND_MAX < rand()) {
	//	//fprintf(stderr," rejected\n", );
	//	//free(ADEnergy_Chaint);
	//	return 0;	/* disregard rejected changes */
	//}

	if (sim_params->protein_model.external_potential_type == 5 && externalloss < 0.0 && !sim_params->NS &&   externalloss * RAND_MAX * external_k < -rand()) {
	//free(ADEnergy_Chaint);
		if (sim_params->protein_model.external_potential_type == 5)
			free(ADEnergy_Chaint);
		return 0;
	}



	
	//if (sim_params->protein_model.external_potential_type == 5 && !sim_params->NS && externalloss < 0.0 &&   exp(sim_params->thermobeta * externalloss *external_k) * RAND_MAX < rand()) {
	//	//fprintf(stderr," rejected\n");
	//	//free(ADEnergy_Chaint);
	//	return 0;	/* disregard rejected changes */
	//	
	//}


	//fprintf(stderr," accepted\n");
	/*Nested Sampling criteria -- important second possibility is for FLEX and is otherwise ignored */



	if(sim_params->NS && ((-logLstar > *currE && -logLstar < *currE - loss) || (-logLstar < *currE && loss < 0  )  )) {
		//free(ADEnergy_Chaint);
		if (sim_params->protein_model.external_potential_type == 5)
			free(ADEnergy_Chaint);
		return 0;
		
	}
	

	/* commit accepted changes */
	for (i = start; i <= end; i++)
		for (j = 1; j < chain->NAA; j++)
	    	chain->Erg(reModNum(i, chain->NAA-1), j) = chain->Erg(j, reModNum(i, chain->NAA-1)) = chaint->Ergt(i, j);
    if (sim_params->protein_model.external_potential_type == 5) {
		for (j = start; j <= end; j++)
			chain->Erg(0, reModNum(j, chain->NAA-1)) = ADEnergy_Chaint[j-start];
		free(ADEnergy_Chaint);
    }
	chain->Erg(0, 0) = 0.0;

	for (j = 1; j < chain->NAA; j++) {
		chain->Erg(0, 0) += chain->Erg(0, j);
	}

	if (sim_params->protein_model.external_potential_type2 == 4) {
		chain->Erg(1, 0) = cyclicBondEnergy;
	}
	*currE -= internalloss + externalloss;
	//free(ADEnergy_Chaint);
    return 1;
}
/* Make a crankshaft move.  This is a local move that involves
the crankshaft rotation of up to 4 peptde bonds.  Propose a
move, and apply the Metropolis criteria. */
void transmutate(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params)
{
	if (sim_params->protein_model.external_potential_type != 5) {
		return 0;
	}
	/*translational move*/
	if (transPtsCount==0) return;
	double transvec[3];
	for (int i = 1; i < chain->NAA; i++) {
		chaint->aat[i] = chain->aa[i];
	}
	casttriplet(chaint->xaat[0], chain->xaa[0]);
	for (int i = 1; i <= chain->NAA - 1; i++) {
		casttriplet(chaint->xaat[i], chain->xaa[i]);
	}
	casttriplet(chaint->xaat_prev[chain->aa[1].chainid], chain->xaa_prev[chain->aa[1].chainid]);

	int transPtsID = rand() % transPtsCount;


	int centerAAID = (chain->NAA + 1) / 2;
	transvec[0] = - chain->aa[centerAAID].c[0] + Xpts[transPtsID];
	transvec[1] = - chain->aa[centerAAID].c[1] + Ypts[transPtsID];
	transvec[2] =  -chain->aa[centerAAID].c[2] + Zpts[transPtsID];
	double movement = 0;
	//fprintf(stderr, "translational move %g %g %g %d \n", transvec[0][vecind], transvec[1][vecind], transvec[2][vecind],chain->NAA);
	for (int i = 0; i < 3; i++) {
		movement = transvec[i];
		for (int j = 1; j < chain->NAA; j++) {
			if (chaint->aat[j].etc & G__){
				chaint->aat[j].g[i] += movement;
			}
			if (chaint->aat[j].etc & G2_){
				chaint->aat[j].g2[i] += movement;
			}

			if (chaint->aat[j].id != 'P') {
				chaint->aat[j].h[i] += movement;
			}
			chaint->aat[j].n[i] += movement;
			chaint->aat[j].ca[i] += movement;
			chaint->aat[j].c[i] += movement;
			chaint->aat[j].o[i] += movement;
			if (chaint->aat[j].id != 'G') {
				chaint->aat[j].cb[i] += movement;
			}

		}
	}


	double* ADEnergy_Chaint;
	//for (int i = 1; i <= chain->NAA - 1; i++){
	//	ADEnergy_Chaint[i - 1] = ADenergy(chaint->aat + i, &(sim_params->protein_model));
	//}

	ADEnergy_Chaint = ADenergyNoClash(1, chain->NAA-1,chain,chaint,&(sim_params->protein_model), 0);

	//for (int i = 1; i <= chain->NAA-1; i++){
	//	//ADEnergy_Chaint[i-1] = chaint->Ergt(0, i);;
	//	externalloss += chain->Erg(0, i) - ADEnergy_Chaint[i-1];
	//}	



	//externalloss += chain->Erg(0, i) - ADEnergy_Chaint[i-start];


	//double transExtEne = global_energy(1, chain->NAA - 1, chain, chaint, biasmap, &(sim_params->protein_model));


	//if loss is negative, it's worse, bad
	//double externalloss = chain->Erg(0, 0) - transExtEne;
	casttriplet(chain->xaa[0], chaint->xaat[0]);

	for (int i = 1; i <= chain->NAA - 1; i++) {
		casttriplet(chain->xaa[i], chaint->xaat[i]);
	}


	chain->Erg(0, 0) = 0.0;
	//if (transExtEne < -30) fprintf(stderr, "committing moved !!\n");
	for (int j = 1; j < chain->NAA; j++) {
	    chain->Erg(0, j) = ADEnergy_Chaint[j - 1];
		chain->Erg(0, 0) += chain->Erg(0, j);
	}

	for (int i = 1; i <= chain->NAA - 1; i++) {
		chain->aa[i] = chaint->aat[i];
	}
	fprintf(stderr, "transmutate!!! %g %g %g\n", chain->aa[centerAAID].c[0], chain->aa[centerAAID].c[1], chain->aa[centerAAID].c[2]);
	fprintf(stderr, "transmutate!!! %g %g %g\n", Xpts[transPtsID], Ypts[transPtsID], Zpts[transPtsID]);
	//copybetween(chain, chaint);
	free(ADEnergy_Chaint);
}


/* Make a translation move. */
static int transmove(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params)
{
	/*translational move*/
	//
	//double transvec[3][chain->NAA - 1];

	if (sim_params->protein_model.external_potential_type != 5) {
		return 0;
	}
	double transvec[3];
	int i;
	for (i = 1; i < chain->NAA; i++) {
		chaint->aat[i] = chain->aa[i];
	}
	//casttriplet(chaint->xaat[0], chain->xaa[0]);
	//for (int i = 1; i <= chain->NAA - 1; i++) {
	//	casttriplet(chaint->xaat[i], chain->xaa[i]);
	//}
	//casttriplet(chaint->xaat_prev[chain->aa[1].chainid], chain->xaa_prev[chain->aa[1].chainid]);

	double movement[3];
	int vecind1 = 0;
	vecind1 = (int)(((double)rand() / RAND_MAX)*(chain->NAA - 1));
	int vecind2 = vecind1;
	while (vecind2 == vecind1) {
		vecind2 = (int)(((double)rand() / RAND_MAX)*(chain->NAA - 1));
	}
	transvec[0] = chain->aa[vecind1].c[0] - chain->aa[vecind2].c[0];
	transvec[1] = chain->aa[vecind1].c[1] - chain->aa[vecind2].c[1];
	transvec[2] = chain->aa[vecind1].c[2] - chain->aa[vecind2].c[2];
	double length = (double)rand() / RAND_MAX;
	//fprintf(stderr, "translational move %g %g %g %d \n", transvec[0][vecind], transvec[1][vecind], transvec[2][vecind],chain->NAA);
	for (i = 0; i < 3; i++) {	
		movement[i] = 0.0;
		if (chain->Erg(0, 0) > 0) {
			if (length > 0.2) {
				movement[i] = 0.4 * rand() / RAND_MAX - 0.2;
			}
			else {
				movement[i] = transvec[i] * length / abs(vecind2 - vecind1);
			}
		}
		else {
			if (length > 0.2) {
				movement[i] = 0.4 * rand() / RAND_MAX - 0.2;
			}
			else {
				movement[i] = transvec[i] * length / abs(vecind2 - vecind1);
			}
		}	
		for (int j = 1; j < chain->NAA; j++) {
			if (chaint->aat[j].etc & G__){
				chaint->aat[j].g[i] += movement[i];
			}
			if (chaint->aat[j].etc & G2_){
				chaint->aat[j].g2[i] += movement[i];
			}
			if (chaint->aat[j].id != 'P') {
				chaint->aat[j].h[i] += movement[i];
			}
			chaint->aat[j].n[i] += movement[i];
			chaint->aat[j].ca[i] += movement[i];
			chaint->aat[j].c[i] += movement[i];
			chaint->aat[j].o[i] += movement[i];
			if (chaint->aat[j].id != 'G') {
				chaint->aat[j].cb[i] += movement[i];
			}

		}
	}

	double externalloss = 0.0;
	double* ADEnergy_Chaint;
	if (sim_params->protein_model.external_potential_type == 5){
		ADEnergy_Chaint = ADenergyNoClash(1, chain->NAA-1,chain,chaint,&(sim_params->protein_model), 0);
		for (i = 1; i <= chain->NAA-1; i++){
			externalloss += chain->Erg(0, i) - ADEnergy_Chaint[i-1];			
		}	
	}
	//if loss is negative, it's worse, bad

	double external_k = 1.0;
	if (sim_params->protein_model.external_potential_type == 5 || sim_params->protein_model.external_potential_type2 == 4)	external_k = sim_params->protein_model.external_k[0];
	//if (externalloss < -10) external_k = 0.09 * external_k;

	if (chain->Erg(0, 0) > 20 || externalloss < -10) external_k = 0.05 * external_k;
	//if (chain->Erg(0, 0) > 50) external_k = 0.2 * external_k;

	//if (moved && allowed(chain, chaint, biasmap, 1, chain->NAA - 1, logLstar, currE, sim_params)) {
	if (externalloss < 0.0 && externalloss * RAND_MAX * external_k < -rand()) {
		free(ADEnergy_Chaint);
		return 0;
	}

	else {

		//casttriplet(chain->xaa_prev[chain->aa[1].chainid], chaint->xaat_prev[chain->aa[1].chainid]);
		//
		//for (int i = 0; i <= chain->NAA - 1; i++) {
		//	casttriplet(chain->xaa[i], chaint->xaat[i]);
		//}


	    chain->Erg(0, 0) = 0.0;

	    for (int j = 1; j < chain->NAA; j++) {
	        chain->Erg(0, j) = ADEnergy_Chaint[j - 1];
	    	chain->Erg(0, 0) += chain->Erg(0, j);
	    }
		
		for (int i = 1; i <= chain->NAA - 1; i++) {
			chain->aa[i] = chaint->aat[i];
		}
		//copybetween(chain, chaint);
		free(ADEnergy_Chaint);
		return 1;
	}
}

/* Make a translation move. */
int transopt(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params)
{
	if (sim_params->protein_model.external_potential_type != 5) {
		return 0;
	}
	/*translational optimization*/
	//
	//double transvec[3][chain->NAA - 1];
	double transvec[3];
	int i;
	for (i = 1; i < chain->NAA; i++) {
		chaint->aat[i] = chain->aa[i];
	}
	//casttriplet(chaint->xaat[0], chain->xaa[0]);
	//for (i = 1; i <= chain->NAA - 1; i++) {
	//	casttriplet(chaint->xaat[i], chain->xaa[i]);
	//}
	//casttriplet(chaint->xaat_prev[chain->aa[1].chainid], chain->xaa_prev[chain->aa[1].chainid]);
	
	double movement[3];
	int noImprovStep = 0;
	double currExtE = 0.0;
	double *currADEnergy, *ADEnergy_Chaint;
	double extE = chain->Erg(0,0);
	int step = 0;
	for (step = 0; step < 30; step++) {		
		
		if (noImprovStep == 5) break;
		for (i = 0; i < 3; i++) {
			if (step == 0) movement[i] = 0.0;
			if (step == 1) movement[i] = 0.4 * rand()/RAND_MAX - 0.2;
			for (int j = 1; j < chain->NAA; j++) {
				if (chaint->aat[j].etc & G__){
					chaint->aat[j].g[i] += movement[i];
				}
				if (chaint->aat[j].etc & G2_){
					chaint->aat[j].g2[i] += movement[i];
				}
				if (chaint->aat[j].id != 'P') {
					chaint->aat[j].h[i] += movement[i];
				}
				chaint->aat[j].n[i] += movement[i];
				chaint->aat[j].ca[i] += movement[i];
				chaint->aat[j].c[i] += movement[i];
				chaint->aat[j].o[i] += movement[i];
				if (chaint->aat[j].id != 'G') {
					chaint->aat[j].cb[i] += movement[i];
				}
			}
		}
		currADEnergy = ADenergyNoClash(1, chain->NAA-1,chain,chaint,&(sim_params->protein_model), 1);
		currExtE = 0.0;
		for (i = 1; i <= chain->NAA-1; i++){
			currExtE += currADEnergy[i-1];
			//fprintf(stderr,"%g !!!", currADEnergy[i-1]);
		}
		//fprintf(stderr,"translation %g \n!!!", currExtE);
		if (currExtE < extE) {
		//find better energy
			extE = currExtE;
			free(ADEnergy_Chaint);
			ADEnergy_Chaint = currADEnergy;
			noImprovStep = 0;
			for (i = 0; i < 3; i++) movement[i] = 0.4 * rand()/RAND_MAX - 0.2;
		} else {
			free(currADEnergy);
		//redo the change and make the step smaller
			for (i = 0; i < 3; i++) {
				for (int j = 1; j < chain->NAA; j++) {
					if (chaint->aat[j].etc & G__){
						chaint->aat[j].g[i] -= movement[i];
					}
					if (chaint->aat[j].etc & G2_){
						chaint->aat[j].g2[i] -= movement[i];
					}
					if (chaint->aat[j].id != 'P') {
						chaint->aat[j].h[i] -= movement[i];
					}
					chaint->aat[j].n[i] -= movement[i];
					chaint->aat[j].ca[i] -= movement[i];
					chaint->aat[j].c[i] -= movement[i];
					chaint->aat[j].o[i] -= movement[i];
					if (chaint->aat[j].id != 'G') {
						chaint->aat[j].cb[i] -= movement[i];
					}
				}
				if (rand()%100 > 20)
					movement[i] = -movement[i]/2;	
			}
			noImprovStep += 1;
		}
		//fprintf(stderr,"translation %g %g %g %g %g %d!!!\n", currExtE, extE ,movement[0] ,movement[1], movement[2], noImprovStep);
	}
	if (extE - chain->Erg(0,0) > -0.00001) return 0;
	//casttriplet(chain->xaa_prev[chain->aa[1].chainid], chaint->xaat_prev[chain->aa[1].chainid]);
	//
	//for (int i = 0; i <= chain->NAA - 1; i++) {
	//	casttriplet(chain->xaa[i], chaint->xaat[i]);
	//}


	chain->Erg(0, 0) = 0.0;

	for (int j = 1; j < chain->NAA; j++) {
	    chain->Erg(0, j) = ADEnergy_Chaint[j - 1];
	 	chain->Erg(0, 0) += chain->Erg(0, j);
	}

	free(ADEnergy_Chaint);

	
	for (int i = 1; i <= chain->NAA - 1; i++) {
		chain->aa[i] = chaint->aat[i];
	}
	//copybetween(chain, chaint);
	return step;
	
}


/* Make a crankshaft move.  This is a local move that involves
   the crankshaft rotation of up to 4 peptde bonds.  Propose a
   move, and apply the Metropolis criteria. */
static int crankshaft(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params)
{	
	int start, end, len, toss;
	double alpha;
	vector a;
	matrix t;
	const double discrete = 2.0 / RAND_MAX;
    
	//if(sim_params->NS){ 
	for (int i = 1; i < chain->NAA; i++){
		chaint->aat[i].etc = chain->aa[i].etc;
		chaint->aat[i].num = chain->aa[i].num;
		chaint->aat[i].id = chain->aa[i].id;
		chaint->aat[i].chainid = chain->aa[i].chainid;
		chaint->aat[i].SCRot = chain->aa[i].SCRot;
	}
	//}

    
   
	/* setup sidechain dihedral angles */
	/* They change with P = 1/4 (unless fixed) */
	if ((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) {
	    if (!(sim_params->protein_model).fix_chi_angles && rand()/(double)RAND_MAX < 0.25) { /* change chi angles */
		for (int i = 1; i < chain->NAA; i++){ 
		    //fprintf(stderr,"chi angles of amino acid %d",i);
		    //fprintf(stderr," %c",chain->aa[i].id);
		    //fprintf(stderr," %g",chain->aa[i].chi1); //would fail for G,A
		    //fprintf(stderr," %g\n",chain->aa[i].chi2); //would fail for all but V,I,T
		    if(chain->aa[i].id != 'G' && chain->aa[i].id != 'A' && chain->aa[i].chi1 != DBL_MAX) {
			chaint->aat[i].chi1 = sidechain_dihedral(chain->aa[i].id, sim_params->protein_model.sidechain_properties);//aa[i].chi1;
		    }
		    if((chain->aa[i].id == 'V' || chain->aa[i].id == 'I' || chain->aa[i].id == 'T') && chain->aa[i].chi2 != DBL_MAX) {
			chaint->aat[i].chi2 = sidechain_dihedral2(chain->aa[i].id,chaint->aat[i].chi1, sim_params->protein_model.sidechain_properties);//aa[i].chi2;
		    }
		}
	    } else {
		for (int i = 1; i < chain->NAA; i++){
			//TODO: we might need this here: if(chain->aa[i].id != 'G' && chain->aa[i].id != 'A') {
			    chaint->aat[i].chi1 = chain->aa[i].chi1;
			    chaint->aat[i].chi2 = chain->aa[i].chi2;
			//}
		} 
	    }
	}

	// Calculate the look-up table of allowed MC moves on the 1st call.
	// This will avoid moves involving residues on more than 1 chain
	//            TODO:                     and fixed atoms.
	if (sim_params->MC_lookup_table == NULL) {
		fprintf(stderr,"creating MC move lookup table.\n");
		// Use the sequence for this, where chain breaks are marked.
		if (sim_params->seq == NULL || sim_params->sequence == NULL) {
			fprintf(stderr,"sim_params->seq: %s\n",sim_params->seq);
			fprintf(stderr,"sim_params->sequence: %s\n",sim_params->sequence);
			stop("sequence is not present in sim_params for MC lookup table calculation\n");
		}
		//allocate memory
        //fprintf(stderr,"allocating memory: 4 * (%d-1+%d) integers.\n",sim_params->NAA,sim_params->Nchains);
		int N = (sim_params->NAA - 1 + sim_params->Nchains);
		sim_params->MC_lookup_table = (int *)malloc(sizeof(int) * 4 * N);
		sim_params->MC_lookup_table_n = (int *)malloc(sizeof(int) * 4);
		if (!sim_params->MC_lookup_table || !sim_params->MC_lookup_table_n) stop("Unable to allocate memory for sim_params->MC_lookup_table.");

		fprintf(stderr,"Sequence:    ");
		for (int i=1; i<sim_params->NAA; i++) {
			fprintf(stderr,"%c",chain->aa[i].id);
		}
		fprintf(stderr,"\n");
		fprintf(stderr,"Fixed:       ");
		for (int i=1; i<sim_params->NAA; i++) {
			if (chain->aa[i].etc & FIXED)
				fprintf(stderr,"x");
			else
				fprintf(stderr," ");
		}
		fprintf(stderr,"\n");
		fprintf(stderr,"Constrained: ");
		for (int i=1; i<sim_params->NAA; i++) {
			if (chain->aa[i].etc & CONSTRAINED)
				fprintf(stderr,"x");
			else
				fprintf(stderr," ");
		}
		fprintf(stderr,"\n");
		fprintf(stderr,"Chain:       ");
		for (int i=1; i<sim_params->NAA; i++) {
			fprintf(stderr,"%d",chain->aa[i].chainid % 10);
		}
		fprintf(stderr,"\n");

		//fprintf(stderr,"lookup table:\n");
		//fill in lookup table
		for (int i=0; i<4; i++) { // loop length can be 0-4
			int next = 0;
			fprintf(stderr,"len %d bonds, Nchains %d:", i+1, sim_params->Nchains);
			int fixed_moves = 0; //counter for moves disallowed due to fixed atoms
			for (int j=1; j<(sim_params->NAA - i); j++){ //start of loop

				//Check if any of the atoms would be fixed
				int any_fixed = 0;
				for (int k=j; k<j+i+1; k++) {
					if ((chain->aa[k].etc & FIXED) && (chain->aa[j].chainid == chain->aa[k].chainid)) any_fixed = 1;
				}
				if (any_fixed) {					
					if (chain->aa[j].chainid == chain->aa[j + i].chainid) {
						fprintf(stderr, "fixed amino acid in %d-%d, skipping", j, j + i + 1);
						
						fixed_moves++;
					}
					//also count the extra move at the beginning of the chain
					if (j == 1) {
						if (chain->aa[j].chainid == chain->aa[j+i].chainid) fixed_moves ++;
					} else {
						if (chain->aa[j].chainid != chain->aa[j-1].chainid && chain->aa[j].chainid == chain->aa[j+i].chainid) fixed_moves ++;
					}
					//also count the extra move at the end of the mid-chain
					if (chain->aa[j].chainid != chain->aa[j + 1].chainid && !chain->aa[j + 1].etc & FIXED) fixed_moves++;
					continue;
				}

				//no fixed atoms, add the move(s).

				//if it's the beginning of the chain, and the chain is long enough, also add the previous amino acid
				if (j == 1) { //first beginning
					if (chain->aa[j].chainid == chain->aa[j+i].chainid) {
						sim_params->MC_lookup_table[i*N+next] = j-1;
						fprintf(stderr,"*%d ",j-1);
						next ++;
					}
				} else { // j > 1, also check chainID differs from previous aa
					if (chain->aa[j].chainid != chain->aa[j-1].chainid && chain->aa[j].chainid == chain->aa[j+i].chainid && !chain->aa[j - i].etc & FIXED) {
						sim_params->MC_lookup_table[i*N+next] = j-1;
						fprintf(stderr,"*%d ",j-1);
						next ++;
					}
				}
				//if it's inside the chain or at the end
				if (chain->aa[j].chainid == chain->aa[j+i].chainid) {
					sim_params->MC_lookup_table[i*N+next] = j;
					fprintf(stderr,"x%d ",j);
					next ++;
				}
			}
			//fprintf(stderr,"(next=%d)",next);
			fprintf(stderr, "%d+%d != (%d-1+(1-%d)*%d\n", next, fixed_moves, sim_params->NAA, i, sim_params->Nchains);
			int N_len = (sim_params->NAA - 1) + (1 - i) * sim_params->Nchains; //number of possibilities for this len
			if ((next + fixed_moves) != N_len) {
				//fprintf(stderr,"%d+%d != (%d-1+(1-%d)*%d\n", next, fixed_moves, sim_params->NAA,i,sim_params->Nchains);
				stop("Something has gone wrong.  Maybe too short chains?\n");
			//} else {
			//	fprintf(stderr,"%d == (%d-1+(1-%d)*%d\n", (next), sim_params->NAA,i,sim_params->Nchains);
			}
			sim_params->MC_lookup_table_n[i] = next; //the number of valid moves
			for (int j = next; j<(sim_params->NAA - 1 + sim_params->Nchains); j++) {
				sim_params->MC_lookup_table[i*N+next] = -1;
				fprintf(stderr,"%d ",-1);
				next ++;
			}
			fprintf(stderr,"\n");
		}
	}

	//stop("Everything is OK.");

//TODO multi-chain protein
	int pivot_around_end = 0;
	int pivot_around_start = 0;

	toss = rand();
	/* segment length */
	len = toss & 0x3;	/* segment length minus one */
	if (len > chain->NAA - 2)
		len = chain->NAA - 2;

	//fprintf(stderr,"MC move len = %d,",len);
	/* amino acids are numbered from 1 to NAA-1 */
	/* segment could start 1 before the first amino acid and end 1 after the last amino acid (pivot moves) or within the chain (crankshaft moves) */

	int N = (sim_params->NAA - 1 + sim_params->Nchains);
	//int N_len = (chain->NAA - 1) + (1 - len) * sim_params->Nchains; //number of possibilities for this len (this has been checked when building the lookup table)
	int N_len = sim_params->MC_lookup_table_n[len]; //the number of valid moves (taking into account fixed amino acids)
	//fprintf(stderr," random number in [ 0, %d ],",N_len-1);



	/* segment start */
	start = sim_params->MC_lookup_table[len*N + (toss >> 2) % N_len ];
	/* segment end */
	if ((sim_params->protein_model).fix_CA_atoms) {
	    end = start + 1;
	} else {
	    end = start + len + 1;
	}
	if (start < 0 || end < 0) { //hit a -1 in the table!
		stop("Something has gone wrong when selecting aminio acids for the MC move.\n");
	}
	//fprintf(stderr,"move residues %d -- %d (len: %d bonds),",start,end,len+1);
	for (int ai=start; ai<end; ai++) {
		if (chain->aa[ai].etc & FIXED) {
			fprintf(stderr,"residues %d -- %d (len: %d bonds),",start,end,len+1);
			fprintf(stderr,"\n%d is fixed\n",ai);
			stop("crankshaft: tried to move fixed amino acid.\n");
		}
	}
	int swappp = 0;
	//for cyclic peptide Gary Hack
	while (sim_params->protein_model.external_potential_type2 == 4 && (start == 0 || end == chain->NAA) && ((toss%50) > chain->Erg(1, 0))) {
		toss = rand();
		len = toss & 0x3;	/* segment length minus one */
		if (len > chain->NAA - 2)
			len = chain->NAA - 2;
		N_len = sim_params->MC_lookup_table_n[len];
		start = sim_params->MC_lookup_table[len*N + (toss >> 2) % N_len];
		if ((sim_params->protein_model).fix_CA_atoms) {
			end = start + 1;
		}
		else {
			end = start + len + 1;
		}
		swappp = 1;
		//fprintf(stderr, "s %d e %d \n", start,end);
	}
	//if (swappp == 1) fprintf(stderr, "s %d e %d \n", start, end);
	//fprintf(stderr, "s2 %d e %d\n", start, end);

	/* pivot or crankshaft */
	if (start == 0) {
		pivot_around_end = 1;
		//fprintf(stderr," pivot around end\n");
		ampl = 1.5 * ampl;
	} else if (end == sim_params->NAA) {
		pivot_around_start = 1;
		ampl = 1.5 * ampl;
		//fprintf(stderr," pivot around start\n");
	} else if (chain->aa[start].chainid != chain->aa[end].chainid) {
		//fprintf(stderr,"  chainid[start] = %d, chainid[end] = %d\n",chain->aa[start].chainid,chain->aa[end].chainid);
		if (len == 0) {
			/* special case for multi-chain protein at chain break for len=0 (2 amino acids) */
			if (rand() & 0x2) {
				pivot_around_start = 1;
				//fprintf(stderr," pivot around start\n");
			} else {
				pivot_around_end = 1;
				//fprintf(stderr," pivot around end\n");
			}
		} else {
			if (chain->aa[start].chainid == chain->aa[start+1].chainid) {
				pivot_around_start = 1;
				//fprintf(stderr," pivot around start\n");
			} else if (chain->aa[end].chainid == chain->aa[end-1].chainid) {
				pivot_around_end = 1;
				//fprintf(stderr," pivot around end\n");
			} else {
				stop("something has gone wrong at the MC move selection\n");
			}
		}
	//} else {// else crankshaft
		//fprintf(stderr,"  chainid[start] = %d, chainid[end] = %d\n",chain->aa[start].chainid,chain->aa[end].chainid);
		//fprintf(stderr," crankshaft\n");
	}

	
	/* setup fixed ends for crankshaft or pivot */
	if (pivot_around_end != 1) { // there is a fixed start site
		casttriplet(chaint->xaat[start], chain->xaa[start]); //TODO: use start - 1 ??
		castvec(chaint->aat[start].ca, chain->aa[start].ca);
		//TODO we will also need xaa[start-1]
		if (start == 1) {//if chain start, use x_prev
			casttriplet(chaint->xaat_prev[chain->aa[end].chainid], chain->xaa_prev[chain->aa[end].chainid]);
			//fprintf(stderr,"copying xaat_prev0 %d \n",chain->aa[end].chainid);
			//fprintf(stderr," %g %g %g %g %g %g %g %g %g\n",chaint->xaat_prev[chain->aa[end].chainid][0][0],chaint->xaat_prev[chain->aa[end].chainid][0][1],chaint->xaat_prev[chain->aa[end].chainid][0][2],chaint->xaat_prev[chain->aa[end].chainid][1][0],chaint->xaat_prev[chain->aa[end].chainid][1][1],chaint->xaat_prev[chain->aa[end].chainid][1][2],chaint->xaat_prev[chain->aa[end].chainid][2][0],chaint->xaat_prev[chain->aa[end].chainid][2][1],chaint->xaat_prev[chain->aa[end].chainid][2][2]);
		} else if (chain->aa[start].chainid != chain->aa[start - 1].chainid) {
			casttriplet(chaint->xaat_prev[chain->aa[end].chainid], chain->xaa_prev[chain->aa[end].chainid]);
			//fprintf(stderr,"copying xaat_prev1 %d\n",chain->aa[end].chainid);
			//fprintf(stderr," %g %g %g %g %g %g %g %g %g\n",chaint->xaat_prev[chain->aa[end].chainid][0][0],chaint->xaat_prev[chain->aa[end].chainid][0][1],chaint->xaat_prev[chain->aa[end].chainid][0][2],chaint->xaat_prev[chain->aa[end].chainid][1][0],chaint->xaat_prev[chain->aa[end].chainid][1][1],chaint->xaat_prev[chain->aa[end].chainid][1][2],chaint->xaat_prev[chain->aa[end].chainid][2][0],chaint->xaat_prev[chain->aa[end].chainid][2][1],chaint->xaat_prev[chain->aa[end].chainid][2][2]);
		} else {
			casttriplet(chaint->xaat[start-1], chain->xaa[start-1]);
		}
	} else {
		//we will also need the xaa[start-1], stored in xaa_prev for chain beginnings
		casttriplet(chaint->xaat_prev[chain->aa[end].chainid], chain->xaa_prev[chain->aa[end].chainid]);
		//fprintf(stderr,"copying xaat_prev2\n");
		//fprintf(stderr," %g %g %g %g %g %g %g %g %g\n",chaint->xaat_prev[chain->aa[end].chainid][0][0],chaint->xaat_prev[chain->aa[end].chainid][0][1],chaint->xaat_prev[chain->aa[end].chainid][0][2],chaint->xaat_prev[chain->aa[end].chainid][1][0],chaint->xaat_prev[chain->aa[end].chainid][1][1],chaint->xaat_prev[chain->aa[end].chainid][1][2],chaint->xaat_prev[chain->aa[end].chainid][2][0],chaint->xaat_prev[chain->aa[end].chainid][2][1],chaint->xaat_prev[chain->aa[end].chainid][2][2]);
	}
	if (pivot_around_start != 1) { // there is a fixed end site
		casttriplet(chaint->xaat[end], chain->xaa[end]);
		castvec(chaint->aat[end].ca, chain->aa[end].ca);
	}



	/* magnitude of rotation */
	/* rotate triplets, alpha in [-ampl; +ampl] */
	alpha = ampl * (discrete * rand() - 1.0);

	/* axis of rotation */
	if (pivot_around_start != 1 && pivot_around_end != 1) {
		/* CA_start->CA_end vector for internal crankshaft */
		subtract(a, chain->aa[end].ca, chain->aa[start].ca);
		normalize(a);
	} else 
		/* random vector for pivot at chain end */
		randvector(a);

	/* rotation matrix */
	rotmatrix(t, a, alpha);


	/* rotating the CA_i->CA_i+1 vectors */
	for (int i = start; i < end; i++){
		if (pivot_around_end == 1 && i == start) {
			//do not change the xaa of the previous chain, use this chain's xaa_prev instead
			rotation(chaint->xaat_prev[chain->aa[end].chainid], t, chain->xaa_prev[chain->aa[end].chainid]);
		} else {
			rotation(chaint->xaat[i], t, chain->xaa[i]);
		}
	}
	/* build trial amino acid CAs using the CA-CA vectors */
	if (pivot_around_end != 1) { // start rotation from the start site
		for (int i = start; i < end - 1; i++){ //moving residues start+1 to end-1
			carbonate_f(chaint->aat + i + 1, chaint->aat + i, chaint->xaat[i]);
		}
		if (pivot_around_start == 1) end --;
	}
	else { //pivot around end
		for (int i = end - 1; i > start; i--){ //moving residues end-1 to start+1
			carbonate_b(chaint->aat + i, chaint->aat + i + 1, chaint->xaat[i]);
		}
		start ++;
	}
	

	//if (swappp == 1) fprintf(stderr, "s %d e %d \n", start, end);

	//building the peptide bonds of the amino acids
	//by now start and end have been adjusted if pivoting
	for (int i = start; i <= end; i++){
		//if starting at the beginning of the chain with pivot or crankshaft
		//fprintf(stderr,"metropolis pivot %d", pivot_around_end);
		//fprintf(stderr," start %d", start);
		//fprintf(stderr," current %d", i);
		//if (pivot_around_end == 1 && i == start){
		//	fprintf(stderr,"\n");
		//} else {
		//	fprintf(stderr," chain(current) %d", chain->aa[i].chainid);
		//	fprintf(stderr," chain(prev) %d\n", chain->aa[i-1].chainid);
		//}
		if ((pivot_around_end == 1 && i == start) || (chain->aa[i].chainid != chain->aa[i-1].chainid))  {
			//use this chain's xaa_prev for the the direction of the N-terminal NH
			//fprintf(stderr,"acidate %d with xaat_prev1 %g %g %g %g %g %g %g %g %g\n",i,chaint->xaat_prev[chain->aa[i].chainid][0][0],chaint->xaat_prev[chain->aa[i].chainid][0][1],chaint->xaat_prev[chain->aa[i].chainid][0][2],chaint->xaat_prev[chain->aa[i].chainid][1][0],chaint->xaat_prev[chain->aa[i].chainid][1][1],chaint->xaat_prev[chain->aa[i].chainid][1][2],chaint->xaat_prev[chain->aa[i].chainid][2][0],chaint->xaat_prev[chain->aa[i].chainid][2][1],chaint->xaat_prev[chain->aa[i].chainid][2][2]);
			//fprintf(stderr,"acidate %d with xaat2 %g %g %g %g %g %g %g %g %g\n",i,chaint->xaat[i][0][0],chaint->xaat[i][0][1],chaint->xaat[i][0][2],chaint->xaat[i][1][0],chaint->xaat[i][1][1],chaint->xaat[i][1][2],chaint->xaat[i][2][0],chaint->xaat[i][2][1],chaint->xaat[i][2][2]);
			acidate(chaint->aat + i, chaint->xaat_prev[chain->aa[i].chainid], chaint->xaat[i], sim_params);
		} else {
			//fprintf(stderr,"acidate %d with xaat1 %g %g %g %g %g %g %g %g %g\n",i,chaint->xaat[i-1][0][0],chaint->xaat[i-1][0][1],chaint->xaat[i-1][0][2],chaint->xaat[i-1][1][0],chaint->xaat[i-1][1][1],chaint->xaat[i-1][1][2],chaint->xaat[i-1][2][0],chaint->xaat[i-1][2][1],chaint->xaat[i-1][2][2]);
			//fprintf(stderr,"acidate %d with xaat2 %g %g %g %g %g %g %g %g %g\n",i,chaint->xaat[i][0][0],chaint->xaat[i][0][1],chaint->xaat[i][0][2],chaint->xaat[i][1][0],chaint->xaat[i][1][1],chaint->xaat[i][1][2],chaint->xaat[i][2][0],chaint->xaat[i][2][1],chaint->xaat[i][2][2]);
			acidate(chaint->aat + i, chaint->xaat[i - 1], chaint->xaat[i], sim_params);
		}
	}



	
        /* testing if move is allowed */
	if (!allowed(chain,chaint,biasmap,start, end, logLstar,currE, sim_params))
		return 0;	/* disregard rejected changes */
	//if (swappp == 1) fprintf(stderr, "s %d e %d \n", start, end);
	
	/* commit accepted changes */
	
	//fprintf(stderr,"committing amino acid xaa %d - %d\n",start-1,end);
	if ((pivot_around_end == 1) || (chain->aa[start].chainid != chain->aa[start-1].chainid))  { //update this chain's xaa_prev
		casttriplet(chain->xaa_prev[chain->aa[end].chainid], chaint->xaat_prev[chain->aa[end].chainid]);
		//fprintf(stderr,"saving chain-%d xaat_prev1 %g %g %g %g %g %g %g %g %g\n",chain->aa[end].chainid,chaint->xaat_prev[chain->aa[i].chainid][0][0],chaint->xaat_prev[chain->aa[i].chainid][0][1],chaint->xaat_prev[chain->aa[i].chainid][0][2],chaint->xaat_prev[chain->aa[i].chainid][1][0],chaint->xaat_prev[chain->aa[i].chainid][1][1],chaint->xaat_prev[chain->aa[i].chainid][1][2],chaint->xaat_prev[chain->aa[i].chainid][2][0],chaint->xaat_prev[chain->aa[i].chainid][2][1],chaint->xaat_prev[chain->aa[i].chainid][2][2]);
	} else {
		casttriplet(chain->xaa[start-1], chaint->xaat[start-1]);
	}
	for (int i = start; i <= end; i++){
		casttriplet(chain->xaa[i], chaint->xaat[i]);
	}





	//fprintf(stderr,"committing amino acid aa %d - %d\n",start,end);
	for (int i = start; i <= end; i++) {
		chain->aa[i] = chaint->aat[i];
	}

	return 1;
}




/* Make a crankshaft move.  This is a local move that involves
   the crankshaft rotation of up to 4 peptde bonds.  Propose a
   move, and apply the Metropolis criteria. */
static int crankshaftcyclic(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params)
{	
	int start, end, len, toss;
	double alpha;
	vector a;
	matrix t;
	const double discrete = 2.0 / RAND_MAX;
    
	//if(sim_params->NS){ 
	for (int j = 1; j < chain->NAA; j++){
		chaint->aat[j].etc = chain->aa[j].etc;
		chaint->aat[j].num = chain->aa[j].num;
		chaint->aat[j].id = chain->aa[j].id;
		chaint->aat[j].chainid = chain->aa[j].chainid;
		chaint->aat[j].SCRot = chain->aa[j].SCRot;
		for(int i = 0; i < 3; i++){
			chaint->aat[j].h[i] = chain->aa[j].h[i];	
			chaint->aat[j].n[i] =  chain->aa[j].n[i];		
			chaint->aat[j].ca[i] = chain->aa[j].ca[i];		
			chaint->aat[j].c[i] = chain->aa[j].c[i];		
			chaint->aat[j].o[i] = chain->aa[j].o[i];
			chaint->aat[j].cb[i] = chain->aa[j].cb[i];
			chaint->aat[j].g[i] = chain->aa[j].g[i];
			chaint->aat[j].g2[i] = chain->aa[j].g2[i];
			chaint->xaat[j][i][0] = chain->xaa[j][i][0];
			chaint->xaat[j][i][1] = chain->xaa[j][i][1];
			chaint->xaat[j][i][2] = chain->xaa[j][i][2];
		}
		//chaint->aat[j] = chain->aa[j];
	}
	//}

    
   
	/* setup sidechain dihedral angles */
	/* They change with P = 1/4 (unless fixed) */
	if ((sim_params->protein_model).use_gamma_atoms != NO_GAMMA) {
	    if (!(sim_params->protein_model).fix_chi_angles && rand()/(double)RAND_MAX < 0.25) { /* change chi angles */
		for (int i = 1; i < chain->NAA; i++){ 
		    //fprintf(stderr,"chi angles of amino acid %d",i);
		    //fprintf(stderr," %c",chain->aa[i].id);
		    //fprintf(stderr," %g",chain->aa[i].chi1); //would fail for G,A
		    //fprintf(stderr," %g\n",chain->aa[i].chi2); //would fail for all but V,I,T
		    if(chain->aa[i].id != 'G' && chain->aa[i].id != 'A' && chain->aa[i].chi1 != DBL_MAX) {
			chaint->aat[i].chi1 = sidechain_dihedral(chain->aa[i].id, sim_params->protein_model.sidechain_properties);//aa[i].chi1;
		    }
		    if((chain->aa[i].id == 'V' || chain->aa[i].id == 'I' || chain->aa[i].id == 'T') && chain->aa[i].chi2 != DBL_MAX) {
			chaint->aat[i].chi2 = sidechain_dihedral2(chain->aa[i].id,chaint->aat[i].chi1, sim_params->protein_model.sidechain_properties);//aa[i].chi2;
		    }
		}
	    } else {
		for (int i = 1; i < chain->NAA; i++){
			//TODO: we might need this here: if(chain->aa[i].id != 'G' && chain->aa[i].id != 'A') {
			    chaint->aat[i].chi1 = chain->aa[i].chi1;
			    chaint->aat[i].chi2 = chain->aa[i].chi2;
			//}
		} 
	    }
	}

	toss = rand();
	/* segment length */
	len = toss & 0x3;	/* segment length minus one */
	if (len > chain->NAA - 2)
		len = chain->NAA - 2;

	/* segment start */
	start = reModNum(rand(), chain->NAA-1);
	/* segment end */
	end = (start + len + 1) ;
		

	//fprintf(stderr,"committing amino acid xaa %d - %d %d\n",start,end,len);

	/* setup fixed ends for crankshaft or pivot */	
	castvec(chaint->aat[start].ca, chain->aa[start].ca);
	casttriplet(chaint->xaat[reModNum(start-1, chain->NAA-1)], chain->xaa[reModNum(start-1, chain->NAA-1)]);
	castvec(chaint->aat[reModNum(end, chain->NAA-1)].ca, chain->aa[reModNum(end, chain->NAA-1)].ca);
	casttriplet(chaint->xaat[reModNum(end, chain->NAA-1)], chain->xaa[reModNum(end, chain->NAA-1)]);
	
	/* magnitude of rotation */
	/* rotate triplets, alpha in [-ampl; +ampl] */
	//ampl = ampl/5;
	alpha = ampl * (discrete * rand() - 1.0);

	/* axis of rotation */
	/* CA_start->CA_end vector for internal crankshaft */
	subtract(a, chain->aa[reModNum(end, chain->NAA-1)].ca, chain->aa[start].ca);
	normalize(a);


	/* rotation matrix */
	rotmatrix(t, a, alpha);


	/* rotating the CA_i->CA_i+1 vectors */
	for (int i = start; i < end; i++){
		rotation(chaint->xaat[reModNum(i, chain->NAA-1)], t, chain->xaa[reModNum(i, chain->NAA-1)]);
	}

	/* build trial amino acid CAs using the CA-CA vectors */

	for (int i = start; i < end-1; i++){ //moving residues start+1 to end-1
		carbonate_f(chaint->aat + reModNum(i + 1, chain->NAA-1), chaint->aat + reModNum(i, chain->NAA-1), chaint->xaat[reModNum(i, chain->NAA-1)]);
	}

	//building the peptide bonds of the amino acids
	//by now start and end have been adjusted if pivoting
	for (int i = start; i <= end; i++) {
		acidate(chaint->aat + reModNum(i, chain->NAA-1), chaint->xaat[reModNum(i - 1, chain->NAA-1)], chaint->xaat[reModNum(i, chain->NAA-1)], sim_params);
	}
	
        /* testing if move is allowed */
	if (!allowed(chain,chaint,biasmap,start, end, logLstar,currE, sim_params))
		return 0;	/* disregard rejected changes */
	//if (swappp == 1) fprintf(stderr, "s %d e %d \n", start, end);
	
	/* commit accepted changes */
	
	//fprintf(stderr,"committing amino acid xaa %d - %d\n",start-1,end);

	casttriplet(chain->xaa[reModNum(start -1 , chain->NAA-1)], chaint->xaat[reModNum(start -1 , chain->NAA-1)]);
	for (int i = start; i <= end; i++){
		casttriplet(chain->xaa[reModNum(i , chain->NAA-1)], chaint->xaat[reModNum(i , chain->NAA-1)]);
	}





	//fprintf(stderr,"committing amino acid aa %d - %d\n",start,end);
	for (int i = start; i <= end; i++) {
		chain->aa[reModNum(i , chain->NAA-1)] = chaint->aat[reModNum(i , chain->NAA-1)];
	}

	return 1;
}


/* MC move wrapper.  Call crankshaft to make an MC move, and calculate the acceptance rate.
   Possibly adjust "negative" amplitudes towards the desired acceptance rate. */
void move(Chain *chain,Chaint *chaint, Biasmap *biasmap, double logLstar, double *currE, int changeamp, simulation_params *sim_params)
{  /*Changed so amplitude does not depend on history of chain 

	changeamp = 0 for normal use
	changeamp = 1 if we are wanting to use the MC move in calculation of new
				amplitude
	changeamp = -1 if we are reseting accept and reject to start the calculation of 
					a new amplitude 

*/
	//static int score = 0
	static int transaccept = 0, reject = 0;    
	
	if (changeamp == -1) { sim_params->accept_counter = 0; sim_params->reject_counter = 0; transaccept = 0; }
	if (sim_params->protein_model.external_potential_type2 == 4 && chain->Erg(1,0)<0.5) {
		if (crankshaftcyclic(chain,chaint,biasmap,sim_params->amplitude,logLstar,currE, sim_params))
			sim_params->accept_counter++; 
	}
	else if (crankshaft(chain,chaint,biasmap,sim_params->amplitude,logLstar,currE, sim_params)) {	/* accepted */	
		sim_params->accept_counter++; 
		//fprintf(stderr, "crankshaft!\n");
		/*if (changeamp && amplitude < 0.0 && ++score > 16 && amplitude > -M_PI) {
			score = 0;
			amplitude *= 1.1;
		}*/
	}
	else if (sim_params->protein_model.external_potential_type == 5 && rand() % 100 < 0 && transmove(chain, chaint, biasmap, sim_params->amplitude, logLstar, currE, sim_params) ) {	/* accepted */
	//	//sim_params->accept_counter++;
		transaccept++;
	//	//fprintf(stderr, "translation!\n");
	//	/*if (changeamp && amplitude < 0.0 && ++score > 16 && amplitude > -M_PI) {
	//	score = 0;
	//	amplitude *= 1.1;
	//	}*/
	}
	else {		/* rejected */
		sim_params->reject_counter++; 
		/*if (changeamp && amplitude < 0.0 && --score < -32) {
			score = 0;
			amplitude *= 0.9;
		}*/
	}
    // used to be 1024!! gary hack
	if (sim_params->accept_counter + sim_params->reject_counter + transaccept == 100000) {
		sim_params->acceptance = sim_params->accept_counter / 100000.;
		if (sim_params->acceptance<0.01) fprintf(stderr, "low acceptance %g! %d\n", sim_params->acceptance, transaccept);
		if(changeamp){
		  if (sim_params->acceptance_rate_tolerance <= 0) stop("The acceptance rate tolerance must be positive.");
		  if (sim_params->acceptance_rate_tolerance >= 1) stop("The acceptance rate tolerance must be smaller than 1.");
		  if (sim_params->amplitude_changing_factor <= 0) stop("The amplitude changing factor must be positive.");
		  if (sim_params->amplitude_changing_factor >= 1) stop("The amplitude changing factor must be smaller than 1.");
		  if(sim_params->amplitude < 0.0 && sim_params->acceptance < sim_params->acceptance_rate - sim_params->acceptance_rate_tolerance) sim_params->amplitude *= sim_params->amplitude_changing_factor;
		  else if(sim_params->acceptance > sim_params->acceptance_rate + sim_params->acceptance_rate_tolerance) sim_params->amplitude /= sim_params->amplitude_changing_factor;	
		  if(sim_params->amplitude < -M_PI) sim_params->amplitude = -M_PI;
		}
		sim_params->accept_counter = 0;
		sim_params->reject_counter = 0;
		transaccept = 0;
	}
	
}

void finalize(Chain *chain, Chaint *chaint, Biasmap *biasmap){

	freemem_chaint(chaint);
	free(chaint);
	freemem_chain(chain); //free amino acid chain and energy matrix
	free(chain);
	biasmap_finalise(biasmap); //free contact map

}


