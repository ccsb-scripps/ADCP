/*
** This is an implementation of model interactions between two amino acids
** as well within a single amino acid. This is a rather simple force-field.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<float.h>
#include<math.h>

#include"error.h"
#include"params.h"
#include"vector.h"
#include"rotation.h"
#include"aadict.h"
#include"peptide.h"
#include"vdw.h"


/***********************************************************/
/****                     CONSTANTS                     ****/
/***********************************************************/

//const int skip_14_vdw=0;


/***********************************************************/
/****               ENERGY  CONTRIBUTIONS               ****/
/***********************************************************/



/* hard cutoff potential:
   Energy contribution of vdW repulsion
     E = 10.4 (RT) * 0.95 p^2 / dist^2,    if    dist2 > 0.95 p^2
     E =  0.0 (RT)    otherwise
   r1: position of atom1
   r2: position of atom2
   Rmin: minimum separation between the two atoms, typically sum of covalent radii */
inline double vdw_hard_cutoff(vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff) {
        double p = Rmin;
	double dis2, sqp;
	
	dis2 = distance(r1,r2);
	sqp = 0.95 * p * p;	// arbitrary adjustment 

    //return (dis2 < sqp) ? 10.4 * sqp / dis2 : 0.0;
    return (dis2 < sqp) ? clash_energy_at_hard_cutoff * sqp / dis2 : 0.0;
    
}

/* Lennard-Jones potential:
   Energy contribution of vdW interaction
     E = max{ clash_energy_at_hard_cutoff ; 4 * depth * ( (Rmin/|r1-r2|)^12 - (Rmin/|r1-r2|)^6 ) },    if    dist2 < (0.4 * Rmin^2)
     E = depth * ( (Rmin/|r1-r2|)^12 - (Rmin/|r1-r2|)^6 ),    if    dist2 > (Rmin * vdw_rel_cutoff)^2
     E =  0.0 (RT)    otherwise
   r1: position of atom1
   r2: position of atom2
   Rmin: the minimum energy separation of the two atoms
   depth: the energy at the minimum energy separation, that is the depth of the potential well
   vdw_rel_cutoff: the cutoff of the potential is (Rmin * vdw_rel_cutoff)
   energy_shift: the energy at the cutoff, to be subtracted to have continuous a function */
inline double vdw_lj(vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff) {

	double dis2, p6;// p12;
	int clashing = 0;
	dis2 = distance(r1,r2);
	double Rmin2 = Rmin * Rmin;
	// apply cutoff
	//fprintf(stderr,"vdw_rel_cutoff: %g",vdw_rel_cutoff);
	//fprintf(stderr,"Rmin2: %g",Rmin2);
	//fprintf(stderr,"vector atom1: %g",r1[0]);
	//fprintf(stderr," %g",r1[1]);
	//fprintf(stderr," %g",r1[2]);
	//fprintf(stderr,"vector atom2: %g",r2[0]);
	//fprintf(stderr," %g",r2[1]);
	//fprintf(stderr," %g",r2[2]);
	//fprintf(stderr," dis2: %g\n",dis2);
	if(dis2 > (vdw_rel_cutoff * vdw_rel_cutoff * Rmin2)) {
//		fprintf(stderr,"return 0 dis2 %g\n",dis2);
		return 0.0;
//	} else {
//		fprintf(stderr,"return non0 dis2 %g\n",dis2);
	}
	// set a max energy at an inner cutoff
	if(dis2 < 0.4 * Rmin2) {
		// return clash value
		clashing = 1; //return clash_energy_at_hard_cutoff;
		dis2 = 0.4 * Rmin2;
	}

	// calculate energy, with energy shifted to 0 at the cutoff
	p6 = Rmin2 / dis2;
	p6 *= (p6 * p6);

	double ans = depth * p6 * (p6 - 2.0) - energy_shift;
	if (clashing) {
		if (ans > clash_energy_at_hard_cutoff)
			return ans;
		else
			return clash_energy_at_hard_cutoff;
	} else {
		return ans;
	}
}

/* vdW energy contribution within the atoms of 1 amino acid */
/* Only calculate contributions between atoms separated by >=4 bonds */
/* Revised by Csilla, 2011-10-28 */
double clash(AA *a, model_params *mod_params)
{
	double erg = 0.0;
	double rg;
	double depth = 0;
	double erg_tmp = 0.0;

	//Assign function pointer to the vdW model
	double (*vdw_fn) (vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff) = NULL;
	if (mod_params->vdw_potential == HARD_CUTOFF_VDW_POTENTIAL) {
		vdw_fn = vdw_hard_cutoff;
	} else if (mod_params->vdw_potential == LJ_VDW_POTENTIAL) {
		vdw_fn = vdw_lj;
	} else {
		stop("Clash cannot be calculated without a valid vdW potential.");
	}

//	/* CB_ -- O__*/
//	if (a->etc & CB_ && a->etc & O__ && mod_params->ro > 0.0 && mod_params->rcb > 0.0) {
//		erg_tmp = vdw_fn(a->cb, a->o, mod_params->rcb + mod_params->ro,mod_params->vdw_depth_cb_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_o, mod_params->vdw_clash_energy_at_hard_cutoff);
//		if (erg_tmp > 5.) fprintf(stderr,"Clash within residue %d CB_-O__.\n",a->num);
//		erg += erg_tmp;
//	}

	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		if (a->etc & G__) {
		   if ((rg = sidechain_vdw_radius(a->id,1,mod_params->sidechain_properties)) < 0.0) {
		      fprintf(stderr,"Negative vdw radius for amino acid %c\n",a->id);
		      //exit(EXIT_FAILURE);
		   } else {
		      /* G__ -- O__ */
		      if (a->etc & O__ && mod_params->ro > 0.0) {
			depth = sidechain_vdw_depth_sqrt(a->id,1,mod_params->sidechain_properties)*mod_params->vdw_depth_o_sqrt;
			erg_tmp = vdw_fn(a->g, a->o, rg + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
//			if (erg_tmp > 5.) fprintf(stderr,"Clash within residue %d G__-O__.\n",a->num);
			erg += erg_tmp;
		      }
//		      /* G__ -- C__ */
//		      if (a->etc & C__ && mod_params->rc > 0.0) {
//			depth = sidechain_vdw_depth_sqrt(a->id,1,mod_params->sidechain_properties)*mod_params->vdw_depth_c_sqrt;
//			erg_tmp = vdw_fn(a->g, a->c, rg + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
//			if (erg_tmp > 5.) fprintf(stderr,"Clash within residue %d G__-C__.\n",a->num);
//			erg += erg_tmp;
//		      }
		   }
		}
		if (a->etc & G2_) { /* I, V, T */
		   if ((rg = sidechain_vdw_radius(a->id,2,mod_params->sidechain_properties)) < 0.0) {
		      fprintf(stderr,"Negative vdw radius for amino acid %c\n",a->id);
		      //exit(EXIT_FAILURE);
		   } else {
		      /* G2_ -- O__ */
		      if (a->etc & O__ && mod_params->ro > 0.0) {
			depth = sidechain_vdw_depth_sqrt(a->id,2,mod_params->sidechain_properties)*mod_params->vdw_depth_o_sqrt;
			erg_tmp = vdw_fn(a->g2, a->o, rg + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
//			if (erg_tmp > 5.) fprintf(stderr,"Clash within residue %d G2_-O__.\n",a->num);
			erg += erg_tmp;
		      }
//		      /* G2_ -- C__ */
//		      if (a->etc & C__ && mod_params->rc > 0.0) {
//			depth = sidechain_vdw_depth_sqrt(a->id,2,mod_params->sidechain_properties)*mod_params->vdw_depth_c_sqrt;
//			erg_tmp = vdw_fn(a->g2, a->c, rg + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
//			if (erg_tmp > 5.) fprintf(stderr,"Clash within residue %d G2_-C__.\n",a->num);
//			erg += erg_tmp;
//		      }
		   }
		}
	}

	return erg;
}

double HHvDW(AA *a, AA *b) 
{
	double HHradii = 2.6;
	double HHstrength = 0.1;
	if (a->id == 'P' || b->id == 'P') return 0.;
	double distHH = distance(a->h, b->h);
	if (distHH >= 9) return 0.;
	double p6 = (HHradii*HHradii/distHH) * (HHradii*HHradii/distHH) *(HHradii*HHradii/distHH);
        double erg=0.;
	if (a->id != 'P' && b->id != 'P' && distHH<9) {
		erg = 1/sqrt(distHH) + HHstrength * (p6*p6-2*p6) - 0.2665;
		//fprintf(stderr, "close exN HH %g\n", erg);
	}
        return erg;	
}

#ifdef LJ_NEIGHBOUR_HARD
/* Energy contribution of all vdW interaction between 2 neighbouring residues */
/* order of neighbors matter, b follows a in the chain */
/* Only calculate contributions between atoms separated by >=4 bonds */
/* Only hard cutoff interactions */
double exclude_neighbor(AA *a, AA *b, model_params *mod_params)
{
	double erg = 0.0;
	double rg_a, rg2_a, rg_b, rg2_b;
	rg_a = rg2_a = rg_b = rg2_b = 0.0;


	//Assign function pointer to the vdW model
	double (*vdw_fn) (vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff) = NULL;
	if (mod_params->vdw_potential == HARD_CUTOFF_VDW_POTENTIAL) {
		vdw_fn = vdw_hard_cutoff;
	} else if (mod_params->vdw_potential == LJ_VDW_POTENTIAL) {
		vdw_fn = vdw_lj;
	} else {
		stop("Clash cannot be calculated without a valid vdW potential.");
	}


	erg += HHvDW(a, b);

	/* collect G__ and G2_ vdW parameters */
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		if (a->etc & G__) {
		   if ((rg_a = sidechain_vdw_radius(a->id,1,mod_params->sidechain_properties)) > 0.0);
		}
		if (a->etc & G2_) {
		   if ((rg2_a = sidechain_vdw_radius(a->id,2,mod_params->sidechain_properties)) > 0.0);
		}
		if (b->etc & G__) {
		   if ((rg_b = sidechain_vdw_radius(b->id,1,mod_params->sidechain_properties)) > 0.0);
		}
		if (b->etc & G2_) {
		   if ((rg2_b = sidechain_vdw_radius(b->id,2,mod_params->sidechain_properties)) > 0.0);
		}
	}

	/* these should not be counted in, since 1-4 interactions */
//	if (mod_params->rn > 0.0 && mod_params->rn > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->n, b->n, mod_params->rn + mod_params->rn,mod_params->vdw_depth_n_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_n_n, mod_params->vdw_clash_energy_at_hard_cutoff);
//	if (mod_params->rc > 0.0 && mod_params->rc > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->c, b->c, mod_params->rc + mod_params->rc,mod_params->vdw_depth_c_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_c, mod_params->vdw_clash_energy_at_hard_cutoff);
//	if (mod_params->rc > 0.0 && mod_params->rcb > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->c, b->cb, mod_params->rc + mod_params->rcb,mod_params->vdw_depth_cb_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_c, mod_params->vdw_clash_energy_at_hard_cutoff);
//	if (mod_params->rcb > 0.0 && mod_params->rn > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->cb, b->n, mod_params->rcb + mod_params->rn,mod_params->vdw_depth_cb_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_n, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* All distances will be counted in, since the LJ vdW potential is long ranged */
//
	/* (a) C__ -- (b) O__ */
	if (mod_params->rc > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->c, b->o, mod_params->rc + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
//
	/* (a) O__ -- (b) CB_ */
	if (mod_params->ro > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') erg += vdw_fn(a->o, b->cb, mod_params->ro + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) O__ -- (b) C__ */
	if (mod_params->ro > 0.0 && mod_params->rc > 0.0) erg += vdw_fn(a->o, b->c, mod_params->ro + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) O__ -- (b) O__ */
	if (mod_params->ro > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->o, b->o, mod_params->ro + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* (a) CA_ -- (b) CB_ */
	if (mod_params->rca > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') erg += vdw_fn(a->ca, b->cb, mod_params->rca + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CA_ -- (b) C__ */
	if (mod_params->rca > 0.0 && mod_params->rc > 0.0) erg += vdw_fn(a->ca, b->c, mod_params->rca + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CA_ -- (b) O__ */
	if (mod_params->rca > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->ca, b->o, mod_params->rca + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* (a) CB_ -- (b) CA_ */
	if (mod_params->rcb > 0.0 && mod_params->rca > 0.0 && a->id != 'G') erg += vdw_fn(a->cb, b->ca, mod_params->rcb + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CB_ -- (b) CB_ */
	if (mod_params->rcb > 0.0 && mod_params->rcb > 0.0 && a->id != 'G' && b->id != 'G') erg += vdw_fn(a->cb, b->cb, mod_params->rcb + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CB_ -- (b) C__ */
	if (mod_params->rcb > 0.0 && mod_params->rc > 0.0 && a->id != 'G') erg += vdw_fn(a->cb, b->c, mod_params->rcb + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CB_ -- (b) O__ */
		if (mod_params->rcb > 0.0 && mod_params->ro > 0.0 && a->id != 'G') erg += vdw_fn(a->cb, b->o, mod_params->rcb + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* (a) N__ -- (b) CA_ */
	if (mod_params->rn > 0.0 && mod_params->rca > 0.0) erg += vdw_fn(a->n, b->ca, mod_params->rn + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) N__ -- (b) CB_ */
	if (mod_params->rn > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') erg += vdw_fn(a->n, b->cb, mod_params->rn + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) N__ -- (b) C__ */
	if (mod_params->rn > 0.0 && mod_params->rc > 0.0) erg += vdw_fn(a->n, b->c, mod_params->rn + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) N__ -- (b) O__ */
	if (mod_params->rn > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->n, b->o, mod_params->rn + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* vdW contribution of gamma atoms */

	if (mod_params->use_gamma_atoms != NO_GAMMA) {

	    if (rg_b > 0.0) {
		/* (a) C__ -- (b) G__ */
		if (mod_params->rc > 0.0) {
		   erg += vdw_fn(a->c, b->g, mod_params->rc + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) O__ -- (b) G__ */
		if (mod_params->ro > 0.0) {
		   erg += vdw_fn(a->o, b->g, mod_params->ro + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CA_ -- (b) G__ */
		if (mod_params->rca > 0.0) {
		   erg += vdw_fn(a->ca, b->g, mod_params->rca + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CB_ -- (b) G__ */
		if (rg_b != 0.0 && mod_params->rcb > 0.0 && a->id != 'G') {
		   erg += vdw_fn(a->cb, b->g, mod_params->rcb + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) N__ -- (b) G__ */
		if (rg_b != 0.0 && mod_params->rn > 0.0) {
		   erg += vdw_fn(a->n, b->g, mod_params->rn + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }

	    if (rg2_b > 0.0) {
		/* (a) C__ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->rc > 0.0) {
		   erg += vdw_fn(a->c, b->g2, mod_params->rc + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) O__ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->ro > 0.0) {
		   erg += vdw_fn(a->o, b->g2, mod_params->ro + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CA_ -- (b) G__ */
		if (mod_params->rca > 0.0) {
		   erg += vdw_fn(a->ca, b->g2, mod_params->rca + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CB_ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->rcb > 0.0 && a->id != 'G') {
		   erg += vdw_fn(a->cb, b->g2, mod_params->rcb + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) N__ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->rn > 0.0) {
		   erg += vdw_fn(a->n, b->g2, mod_params->rn + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }

	    if (rg_a != 0.0) {
		/* (a) G__ -- (b) CA_ */
		if (mod_params->rca > 0.0) {
		   erg += vdw_fn(a->g, b->ca, rg_a + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) C__ */
		if (mod_params->rc > 0.0) {
		   erg += vdw_fn(a->g, b->c, rg_a + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) O__ */
		if (mod_params->ro > 0.0) {
		   erg += vdw_fn(a->g, b->o, rg_a + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) CB_ */
		if (mod_params->rcb > 0.0 && b->id != 'G') {
		   erg += vdw_fn(a->g, b->cb, rg_a + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) N__ */
		if (mod_params->rn > 0.0) {
		   erg += vdw_fn(a->g, b->n, rg_a + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) G__ */
		if (rg_b != 0.0) {
		   erg += vdw_fn(a->g, b->g, rg_a + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		if (rg2_b != 0.0) {
		   erg += vdw_fn(a->g, b->g2, rg_a + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }

	    if (rg2_a != 0.0) {
		/* (a) G__ -- (b) CA_ */
		if (mod_params->rca > 0.0) {
		   erg += vdw_fn(a->g2, b->ca, rg2_a + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) C__ */
		if (mod_params->rc > 0.0) {
		   erg += vdw_fn(a->g2, b->c, rg2_a + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) O__ */
		if (mod_params->ro > 0.0) {
		   erg += vdw_fn(a->g2, b->o, rg2_a + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) CB_ */
		if (mod_params->rcb > 0.0 && b->id != 'G') {
		   erg += vdw_fn(a->g2, b->cb, rg2_a + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) N__ */
		if (mod_params->rn > 0.0) {
		   erg += vdw_fn(a->g2, b->n, rg2_a + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) G__ */
		if (rg_b != 0.0) {
		   erg += vdw_fn(a->g2, b->g, rg2_a + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		if (rg2_b != 0.0) {
		   erg += vdw_fn(a->g2, b->g2, rg2_a + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }
	}

	return erg;
}
#else
/* Energy contribution of all vdW interaction between 2 neighbouring residues */
/* order of neighbors matter, b follows a in the chain */
/* Only calculate contributions between atoms separated by >=4 bonds */
double exclude_neighbor(AA *a, AA *b, model_params *mod_params)
{
	double erg = 0.0;
	double rg_a, rg2_a, rg_b, rg2_b;
	rg_a = rg2_a = rg_b = rg2_b = 0.0;
	double eps_g_a_sqrt, eps_g2_a_sqrt, eps_g_b_sqrt, eps_g2_b_sqrt;
	eps_g_a_sqrt =  eps_g2_a_sqrt =  eps_g_b_sqrt =  eps_g2_b_sqrt = 0.0;
	double depth = 0;


	//Assign function pointer to the vdW model
	double (*vdw_fn) (vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff) = NULL;
	if (mod_params->vdw_potential == HARD_CUTOFF_VDW_POTENTIAL) {
		vdw_fn = vdw_hard_cutoff;
	} else if (mod_params->vdw_potential == LJ_VDW_POTENTIAL) {
		vdw_fn = vdw_lj;
	} else {
		stop("Clash cannot be calculated without a valid vdW potential.");
	}

	erg += HHvDW(a, b);

	/* collect G__ and G2_ vdW parameters */
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		if (a->etc & G__) {
		   if ((rg_a = sidechain_vdw_radius(a->id,1,mod_params->sidechain_properties)) > 0.0)
		   eps_g_a_sqrt = sidechain_vdw_depth_sqrt(a->id,1,mod_params->sidechain_properties);
		}
		if (a->etc & G2_) {
		   if ((rg2_a = sidechain_vdw_radius(a->id,2,mod_params->sidechain_properties)) > 0.0)
		   eps_g2_a_sqrt = sidechain_vdw_depth_sqrt(a->id,2,mod_params->sidechain_properties);
		}
		if (b->etc & G__) {
		   if ((rg_b = sidechain_vdw_radius(b->id,1,mod_params->sidechain_properties)) > 0.0)
		   eps_g_b_sqrt = sidechain_vdw_depth_sqrt(b->id,1,mod_params->sidechain_properties);
		}
		if (b->etc & G2_) {
		   if ((rg2_b = sidechain_vdw_radius(b->id,2,mod_params->sidechain_properties)) > 0.0)
		   eps_g2_b_sqrt = sidechain_vdw_depth_sqrt(b->id,2,mod_params->sidechain_properties);
		}
	}

	/* these should not be counted in, since 1-4 interactions */
//	if (mod_params->rn > 0.0 && mod_params->rn > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->n, b->n, mod_params->rn + mod_params->rn,mod_params->vdw_depth_n_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_n_n, mod_params->vdw_clash_energy_at_hard_cutoff);
//	if (mod_params->rc > 0.0 && mod_params->rc > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->c, b->c, mod_params->rc + mod_params->rc,mod_params->vdw_depth_c_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_c, mod_params->vdw_clash_energy_at_hard_cutoff);
//	if (mod_params->rc > 0.0 && mod_params->rcb > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->c, b->cb, mod_params->rc + mod_params->rcb,mod_params->vdw_depth_cb_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_c, mod_params->vdw_clash_energy_at_hard_cutoff);
//	if (mod_params->rcb > 0.0 && mod_params->rn > 0.0 && skip_14_vdw==0) erg += vdw_fn(a->cb, b->n, mod_params->rcb + mod_params->rn,mod_params->vdw_depth_cb_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_n, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* All distances will be counted in, since the LJ vdW potential is long ranged */
//
	/* (a) C__ -- (b) O__ */
	if (mod_params->rc > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->c, b->o, mod_params->rc + mod_params->ro,mod_params->vdw_depth_c_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_o, mod_params->vdw_clash_energy_at_hard_cutoff);
//
	/* (a) O__ -- (b) CB_ */
	if (mod_params->ro > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') erg += vdw_fn(a->o, b->cb, mod_params->ro + mod_params->rcb,mod_params->vdw_depth_cb_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_o, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) O__ -- (b) C__ */
	if (mod_params->ro > 0.0 && mod_params->rc > 0.0) erg += vdw_fn(a->o, b->c, mod_params->ro + mod_params->rc,mod_params->vdw_depth_c_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_o, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) O__ -- (b) O__ */
	if (mod_params->ro > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->o, b->o, mod_params->ro + mod_params->ro,mod_params->vdw_depth_o_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_o_o, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* (a) CA_ -- (b) CB_ */
	if (mod_params->rca > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') erg += vdw_fn(a->ca, b->cb, mod_params->rca + mod_params->rcb,mod_params->vdw_depth_ca_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_cb, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CA_ -- (b) C__ */
	if (mod_params->rca > 0.0 && mod_params->rc > 0.0) erg += vdw_fn(a->ca, b->c, mod_params->rca + mod_params->rc,mod_params->vdw_depth_ca_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_c, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CA_ -- (b) O__ */
	if (mod_params->rca > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->ca, b->o, mod_params->rca + mod_params->ro,mod_params->vdw_depth_ca_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_o, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* (a) CB_ -- (b) CA_ */
	if (mod_params->rcb > 0.0 && mod_params->rca > 0.0 && a->id != 'G') erg += vdw_fn(a->cb, b->ca, mod_params->rcb + mod_params->rca,mod_params->vdw_depth_ca_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_cb, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CB_ -- (b) CB_ */
	if (mod_params->rcb > 0.0 && mod_params->rcb > 0.0 && a->id != 'G' && b->id != 'G') erg += vdw_fn(a->cb, b->cb, mod_params->rcb + mod_params->rcb,mod_params->vdw_depth_cb_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_cb, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CB_ -- (b) C__ */
	if (mod_params->rcb > 0.0 && mod_params->rc > 0.0 && a->id != 'G') erg += vdw_fn(a->cb, b->c, mod_params->rcb + mod_params->rc,mod_params->vdw_depth_cb_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_c, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) CB_ -- (b) O__ */
		if (mod_params->rcb > 0.0 && mod_params->ro > 0.0 && a->id != 'G') erg += vdw_fn(a->cb, b->o, mod_params->rcb + mod_params->ro,mod_params->vdw_depth_cb_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_o, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* (a) N__ -- (b) CA_ */
	if (mod_params->rn > 0.0 && mod_params->rca > 0.0) erg += vdw_fn(a->n, b->ca, mod_params->rn + mod_params->rca,mod_params->vdw_depth_ca_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_n, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) N__ -- (b) CB_ */
	if (mod_params->rn > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') erg += vdw_fn(a->n, b->cb, mod_params->rn + mod_params->rcb,mod_params->vdw_depth_cb_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_n, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) N__ -- (b) C__ */
	if (mod_params->rn > 0.0 && mod_params->rc > 0.0) erg += vdw_fn(a->n, b->c, mod_params->rn + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
	/* (a) N__ -- (b) O__ */
	if (mod_params->rn > 0.0 && mod_params->ro > 0.0) erg += vdw_fn(a->n, b->o, mod_params->rn + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);

	/* vdW contribution of gamma atoms */

	if (mod_params->use_gamma_atoms != NO_GAMMA) {

	    if (rg_b > 0.0) {
		/* (a) C__ -- (b) G__ */
		if (mod_params->rc > 0.0) {
		   depth = eps_g_b_sqrt*mod_params->vdw_depth_c_sqrt;
		   erg += vdw_fn(a->c, b->g, mod_params->rc + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) O__ -- (b) G__ */
		if (mod_params->ro > 0.0) {
		   depth = eps_g_b_sqrt*mod_params->vdw_depth_o_sqrt;
		   erg += vdw_fn(a->o, b->g, mod_params->ro + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CA_ -- (b) G__ */
		if (mod_params->rca > 0.0) {
		   depth = eps_g_b_sqrt*mod_params->vdw_depth_ca_sqrt;
		   erg += vdw_fn(a->ca, b->g, mod_params->rca + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CB_ -- (b) G__ */
		if (rg_b != 0.0 && mod_params->rcb > 0.0 && a->id != 'G') {
		   depth = eps_g_b_sqrt*mod_params->vdw_depth_cb_sqrt;
		   erg += vdw_fn(a->cb, b->g, mod_params->rcb + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) N__ -- (b) G__ */
		if (rg_b != 0.0 && mod_params->rn > 0.0) {
		   depth = eps_g_b_sqrt*mod_params->vdw_depth_n_sqrt;
		   erg += vdw_fn(a->n, b->g, mod_params->rn + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }

	    if (rg2_b > 0.0) {
		/* (a) C__ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->rc > 0.0) {
		   depth = eps_g2_b_sqrt*mod_params->vdw_depth_c_sqrt;
		   erg += vdw_fn(a->c, b->g2, mod_params->rc + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) O__ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->ro > 0.0) {
		   depth = eps_g2_b_sqrt*mod_params->vdw_depth_o_sqrt;
		   erg += vdw_fn(a->o, b->g2, mod_params->ro + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CA_ -- (b) G__ */
		if (mod_params->rca > 0.0) {
		   depth = eps_g2_b_sqrt*mod_params->vdw_depth_ca_sqrt;
		   erg += vdw_fn(a->ca, b->g2, mod_params->rca + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) CB_ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->rcb > 0.0 && a->id != 'G') {
		   depth = eps_g2_b_sqrt*mod_params->vdw_depth_cb_sqrt;
		   erg += vdw_fn(a->cb, b->g2, mod_params->rcb + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) N__ -- (b) G__ */
		if (rg2_b != 0.0 && mod_params->rn > 0.0) {
		   depth = eps_g2_b_sqrt*mod_params->vdw_depth_n_sqrt;
		   erg += vdw_fn(a->n, b->g2, mod_params->rn + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }

	    if (rg_a != 0.0) {
		/* (a) G__ -- (b) CA_ */
		if (mod_params->rca > 0.0) {
		   depth = eps_g_a_sqrt*mod_params->vdw_depth_ca_sqrt;
		   erg += vdw_fn(a->g, b->ca, rg_a + mod_params->rca,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) C__ */
		if (mod_params->rc > 0.0) {
		   depth = eps_g_a_sqrt*mod_params->vdw_depth_c_sqrt;
		   erg += vdw_fn(a->g, b->c, rg_a + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) O__ */
		if (mod_params->ro > 0.0) {
		   depth = eps_g_a_sqrt*mod_params->vdw_depth_o_sqrt;
		   erg += vdw_fn(a->g, b->o, rg_a + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) CB_ */
		if (mod_params->rcb > 0.0 && b->id != 'G') {
		   depth = eps_g_a_sqrt*mod_params->vdw_depth_cb_sqrt;
		   erg += vdw_fn(a->g, b->cb, rg_a + mod_params->rcb,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) N__ */
		if (mod_params->rn > 0.0) {
		   depth = eps_g_a_sqrt*mod_params->vdw_depth_n_sqrt;
		   erg += vdw_fn(a->g, b->n, rg_a + mod_params->rn,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) G__ */
		if (rg_b != 0.0) {
		   depth = eps_g_a_sqrt*eps_g_b_sqrt;
		   erg += vdw_fn(a->g, b->g, rg_a + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		if (rg2_b != 0.0) {
		   depth = eps_g_a_sqrt*eps_g2_b_sqrt;
		   erg += vdw_fn(a->g, b->g2, rg_a + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }

	    if (rg2_a != 0.0) {
		/* (a) G__ -- (b) CA_ */
		if (mod_params->rca > 0.0) {
		   depth = eps_g2_a_sqrt*mod_params->vdw_depth_ca_sqrt;
		   erg += vdw_fn(a->g2, b->ca, rg2_a + mod_params->rca,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) C__ */
		if (mod_params->rc > 0.0) {
		   depth = eps_g2_a_sqrt*mod_params->vdw_depth_c_sqrt;
		   erg += vdw_fn(a->g2, b->c, rg2_a + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) O__ */
		if (mod_params->ro > 0.0) {
		   depth = eps_g2_a_sqrt*mod_params->vdw_depth_o_sqrt;
		   erg += vdw_fn(a->g2, b->o, rg2_a + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) CB_ */
		if (mod_params->rcb > 0.0 && b->id != 'G') {
		   depth = eps_g2_a_sqrt*mod_params->vdw_depth_cb_sqrt;
		   erg += vdw_fn(a->g2, b->cb, rg2_a + mod_params->rcb,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) N__ */
		if (mod_params->rn > 0.0) {
		   depth = eps_g2_a_sqrt*mod_params->vdw_depth_n_sqrt;
		   erg += vdw_fn(a->g2, b->n, rg2_a + mod_params->rn,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		/* (a) G__ -- (b) G__ */
		if (rg_b != 0.0) {
		   depth = eps_g2_a_sqrt*eps_g_b_sqrt;
		   erg += vdw_fn(a->g2, b->g, rg2_a + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
		if (rg2_b != 0.0) {
		   depth = eps_g2_a_sqrt*eps_g2_b_sqrt;
		   erg += vdw_fn(a->g2, b->g2, rg2_a + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
		}
	    }
	}

	return erg;
}
#endif

#ifdef LJ_HBONDED_HARD
/* Energy contribution of all vdW interaction between 2 residues
   a, b: amino acids
   d2: squared distance of their alpha carbons
   hbond_proxmity: if they or their neighbours connect these guys by H-bonds */
double exclude_hard(AA *a, AA *b, double d2, model_params *mod_params, int hbond_proximity)
{
	double erg = 0.0;
	double rg_a, rg2_a, rg_b, rg2_b;
	rg_a =  rg2_a =  rg_b =  rg2_b = 0.0;
	double vdw_erg = 0;
	const double backbone_constants[3] = { mod_params->vdw_backbone_cutoff, mod_params->vdw_backbone_cutoff, mod_params->vdw_backbone_cutoff };

	//Assign function pointer to the vdW model
	double (*vdw_fn) (vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff) = NULL;
	if (mod_params->vdw_potential == HARD_CUTOFF_VDW_POTENTIAL) {
		vdw_fn = vdw_hard_cutoff;
	} else if (mod_params->vdw_potential == LJ_VDW_POTENTIAL) {
		vdw_fn = vdw_lj;
	} else {
		stop("Clash cannot be calculated without a valid vdW potential.");
	}

	erg += HHvDW(a, b);

    /* Calculates the correct index for maxvdw_gamma_gamma 
     * note the other d2 > should be changed for maximum efficiency*/
    int index = (a->id - '@') * 26 + (b->id - '@') -1;

	if (mod_params->use_gamma_atoms != NO_GAMMA) {

		/* 1. check maximum interaction distance: this is the G-G case */
		if (d2 > (mod_params->vdw_gamma_gamma_cutoff)[index]) /* cg - cg */
			return erg;
        
			if (a->etc & G__) {
			   if ((rg_a = sidechain_vdw_radius(a->id,1,mod_params->sidechain_properties)) > 0.0);
			}
			if (a->etc & G2_) {
			   if ((rg2_a = sidechain_vdw_radius(a->id,2,mod_params->sidechain_properties)) > 0.0);
			}
			if (b->etc & G__) {
			   if ((rg_b = sidechain_vdw_radius(b->id,1,mod_params->sidechain_properties)) > 0.0);
			}
			if (b->etc & G2_) {
			   if ((rg2_b = sidechain_vdw_radius(b->id,2,mod_params->sidechain_properties)) > 0.0);
			}
			

		if (rg_a != 0.0) {
			/* Skip CYS--CYS gamma-gamma interactions (should really only miss when they are S-S bonded!, but let's try this for now) */
			if (rg_b != 0.0 && (a->id != 'C' || b->id != 'C')) { //&& (mod_params->Sbond_strength == 0 || a->id != 'C' || b->id != 'C' )) 
			  vdw_erg = vdw_fn(a->g, b->g, rg_a + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
			if (rg2_b != 0.0) {
			  vdw_erg = vdw_fn(a->g, b->g2, rg_a + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
		}
		if (rg2_a != 0.0) {
			if (rg_b != 0.0) {
			   vdw_erg = vdw_fn(a->g2, b->g, rg2_a + rg_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
			if (rg2_b != 0.0) {
			   vdw_erg = vdw_fn(a->g2, b->g2, rg2_a + rg2_b,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
		}
        
		/* 2. check slightly closer: this is the G-nonG case */
		if (d2 > (mod_params->vdw_gamma_nongamma_cutoff)[index]) /* cg - o */
			return erg;
        
		if (rg_a != 0.0 && mod_params->ro > 0.0) {
			vdw_erg = vdw_fn(a->g, b->o, rg_a + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->ro > 0.0) {
			vdw_erg = vdw_fn(a->g2, b->o, rg2_a + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->ro > 0.0) {
			vdw_erg = vdw_fn(b->g, a->o, rg_b + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->ro > 0.0) {
			vdw_erg = vdw_fn(b->g2, a->o, rg2_b + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 94.79) /* cg - cb */
		//	return erg;
        
		/* Skip CYS--CYS beta-gamma interactions (should really only miss when they are S-S bonded!, but let's try this for now) */
		if (rg_a != 0.0 && mod_params->rcb > 0.0 && (a->id != 'C' || b->id != 'C') && b->id != 'G') {
			vdw_erg = vdw_fn(a->g, b->cb, rg_a + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
			vdw_erg = vdw_fn(a->g2, b->cb, rg2_a + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rcb > 0.0 && (a->id != 'C' || b->id != 'C') && a->id != 'G') {
			vdw_erg = vdw_fn(b->g, a->cb, rg_b + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rcb > 0.0 && a->id != 'G') {
			vdw_erg = vdw_fn(b->g2, a->cb, rg2_b + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 91.80) /* cg -- c */
		//	return erg;
        
		if (rg_a != 0.0 && mod_params->rc > 0.0) {
			vdw_erg = vdw_fn(a->g, b->c, rg_a + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rc > 0.0) {
			vdw_erg = vdw_fn(a->g2, b->c, rg2_a + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rc > 0.0) {
			vdw_erg = vdw_fn(b->g, a->c, rg_b + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rc > 0.0) {
			vdw_erg = vdw_fn(b->g2, a->c, rg2_b + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 88.06) /* cg -- n */
		//	return erg;
        
		if (rg_a != 0.0 && mod_params->rn > 0.0) {
			vdw_erg = vdw_fn(a->g, b->n, rg_a + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rn > 0.0) {
			vdw_erg = vdw_fn(a->g2, b->n, rg2_a + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rn > 0.0) {
			vdw_erg = vdw_fn(b->g, a->n, rg_b + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rn > 0.0) {
			vdw_erg = vdw_fn(b->g2, a->n, rg2_b + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 67.33) /* cg - ca */
		//	return erg;
        
		if (rg_a != 0.0 && mod_params->rca > 0.0) {
			vdw_erg = vdw_fn(a->g, b->ca, rg_a + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rca > 0.0) {
			vdw_erg = vdw_fn(a->g2, b->ca, rg2_a + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rca > 0.0) {
			vdw_erg = vdw_fn(b->g, a->ca, rg_b + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rca > 0.0) {
			vdw_erg = vdw_fn(b->g2, a->ca, rg2_b + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}

	}

	/* 3. check even closer: this is the nonG-nonG case (CB-extended backbone) */
	if (d2 > backbone_constants[0]) /* o-o, farthest */
	/* 3.A  O-else */
		return erg;
    /* o - o, o - n, o - c */

	if (mod_params->ro > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->o, b->o, mod_params->ro + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	// only hard cutoff for N-O, the H-bond interactions will pull them in
	if (mod_params->ro > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->o, b->n, mod_params->ro + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->rn > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->n, b->o, mod_params->rn + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->ro > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->o, b->c, mod_params->ro + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->c, b->o, mod_params->rc + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}

   if ( d2 > backbone_constants[1])
		return erg; 

	/* 3.B  any-any, within N-N, mostly N with others */
	/* n - n, n - c, c - c, o - cb, o - ca, n - cb, c - cb */
	if (mod_params->rn > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->n, b->n, mod_params->rn + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rn > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->n, b->c, mod_params->rn + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->c, b->n, mod_params->rc + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->c, b->c, mod_params->rc + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->ro > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->o, b->cb, mod_params->ro + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->ro > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->o, mod_params->rcb + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->ro > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->o, b->ca, mod_params->ro + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->o, mod_params->rca + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rn > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->n, b->cb, mod_params->rn + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rn > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->n, mod_params->rcb + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->c, b->cb, mod_params->rc + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rc > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->c, mod_params->rcb + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	
	if ( d2 > backbone_constants[2])
		return erg; 
	
	/* 3.B  rest, within N-CA, mostly ca with all others */
	/* n - ca, c - ca, cb - ca, cb - cb, ca - ca */
	if (mod_params->rca > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->n, mod_params->rca + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rn > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->n, b->ca, mod_params->rn + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->c, b->ca, mod_params->rc + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->c, mod_params->rca + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rcb > 0.0 && a->id != 'G' && b->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->cb, mod_params->rcb + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->ca, b->cb, mod_params->rca + mod_params->rcb,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rca > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->ca, mod_params->rcb + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->ca, mod_params->rca + mod_params->rca,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
    
	return erg;
}
#endif
/* Energy contribution of all vdW interaction between 2 residues
   a, b: amino acids
   d2: squared distance of their alpha carbons */
double exclude(AA *a, AA *b, double d2, model_params *mod_params)
{
	double erg = 0.0;
	double rg_a, rg2_a, rg_b, rg2_b;
	rg_a =  rg2_a =  rg_b =  rg2_b = 0.0;
	double eps_g_a_sqrt, eps_g2_a_sqrt, eps_g_b_sqrt, eps_g2_b_sqrt;
	eps_g_a_sqrt =  eps_g2_a_sqrt =  eps_g_b_sqrt =  eps_g2_b_sqrt = 0.0;
	double depth = 0;
	double vdw_erg = 0;
	const double backbone_constants[3] = { mod_params->vdw_backbone_cutoff, mod_params->vdw_backbone_cutoff, mod_params->vdw_backbone_cutoff };

	erg += HHvDW(a, b);

	//Assign function pointer to the vdW model
	double (*vdw_fn) (vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_energy_at_hard_cutoff) = NULL;
	if (mod_params->vdw_potential == HARD_CUTOFF_VDW_POTENTIAL) {
		vdw_fn = vdw_hard_cutoff;
	} else if (mod_params->vdw_potential == LJ_VDW_POTENTIAL) {
		vdw_fn = vdw_lj;
	} else {
		stop("Clash cannot be calculated without a valid vdW potential.");
	}

    /* Calculates the correct index for maxvdw_gamma_gamma 
     * note the other d2 > should be changed for maximum efficiency*/
    int index = (a->id - '@') * 26 + (b->id - '@') -1;

	if (mod_params->use_gamma_atoms != NO_GAMMA) {

		/* 1. check maximum interaction distance: this is the G-G case */
		if (d2 > (mod_params->vdw_gamma_gamma_cutoff)[index]) /* cg - cg */
			return erg;
        
			if (a->etc & G__) {
			   if ((rg_a = sidechain_vdw_radius(a->id,1,mod_params->sidechain_properties)) > 0.0)
			   eps_g_a_sqrt = sidechain_vdw_depth_sqrt(a->id,1,mod_params->sidechain_properties);
			}
			if (a->etc & G2_) {
			   if ((rg2_a = sidechain_vdw_radius(a->id,2,mod_params->sidechain_properties)) > 0.0)
			   eps_g2_a_sqrt = sidechain_vdw_depth_sqrt(a->id,2,mod_params->sidechain_properties);
			}
			if (b->etc & G__) {
			   if ((rg_b = sidechain_vdw_radius(b->id,1,mod_params->sidechain_properties)) > 0.0)
			   eps_g_b_sqrt = sidechain_vdw_depth_sqrt(b->id,1,mod_params->sidechain_properties);
			}
			if (b->etc & G2_) {
			   if ((rg2_b = sidechain_vdw_radius(b->id,2,mod_params->sidechain_properties)) > 0.0)
			   eps_g2_b_sqrt = sidechain_vdw_depth_sqrt(b->id,2,mod_params->sidechain_properties);
			}
			

		if (rg_a != 0.0) {
			/* Skip CYS--CYS gamma-gamma interactions (should really only miss when they are S-S bonded!, but let's try this for now) */
			if (rg_b != 0.0 && (a->id != 'C' || b->id != 'C')) { //&& (mod_params->Sbond_strength == 0 || a->id != 'C' || b->id != 'C' )) 
			  depth = eps_g_a_sqrt*eps_g_b_sqrt;
			  vdw_erg = vdw_fn(a->g, b->g, rg_a + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
			if (rg2_b != 0.0) {
			  depth = eps_g_a_sqrt*eps_g2_b_sqrt;
			  vdw_erg = vdw_fn(a->g, b->g2, rg_a + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
		}
		if (rg2_a != 0.0) {
			if (rg_b != 0.0) {
			  depth = eps_g2_a_sqrt*eps_g_b_sqrt;
			   vdw_erg = vdw_fn(a->g2, b->g, rg2_a + rg_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
			if (rg2_b != 0.0) {
			  depth = eps_g2_a_sqrt*eps_g2_b_sqrt;
			   vdw_erg = vdw_fn(a->g2, b->g2, rg2_a + rg2_b,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
			}
		}
        
		/* 2. check slightly closer: this is the G-nonG case */
		if (d2 > (mod_params->vdw_gamma_nongamma_cutoff)[index]) /* cg - o */
			return erg;
        
		if (rg_a != 0.0 && mod_params->ro > 0.0) {
			depth = eps_g_a_sqrt*mod_params->vdw_depth_o_sqrt;
			vdw_erg = vdw_fn(a->g, b->o, rg_a + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->ro > 0.0) {
			depth = eps_g2_a_sqrt*mod_params->vdw_depth_o_sqrt;
			vdw_erg = vdw_fn(a->g2, b->o, rg2_a + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->ro > 0.0) {
			depth = eps_g_b_sqrt*mod_params->vdw_depth_o_sqrt;
			vdw_erg = vdw_fn(b->g, a->o, rg_b + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->ro > 0.0) {
			depth = eps_g2_b_sqrt*mod_params->vdw_depth_o_sqrt;
			vdw_erg = vdw_fn(b->g2, a->o, rg2_b + mod_params->ro,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 94.79) /* cg - cb */
		//	return erg;
        
		/* Skip CYS--CYS beta-gamma interactions (should really only miss when they are S-S bonded!, but let's try this for now) */
		if (rg_a != 0.0 && mod_params->rcb > 0.0 && (a->id != 'C' || b->id != 'C') && b->id != 'G') {
			depth = eps_g_a_sqrt*mod_params->vdw_depth_cb_sqrt;
			vdw_erg = vdw_fn(a->g, b->cb, rg_a + mod_params->rcb,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
			depth = eps_g2_a_sqrt*mod_params->vdw_depth_cb_sqrt;
			vdw_erg = vdw_fn(a->g2, b->cb, rg2_a + mod_params->rcb,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rcb > 0.0 && (a->id != 'C' || b->id != 'C') && a->id != 'G') {
			depth = eps_g_b_sqrt*mod_params->vdw_depth_cb_sqrt;
			vdw_erg = vdw_fn(b->g, a->cb, rg_b + mod_params->rcb,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rcb > 0.0 && a->id != 'G') {
			depth = eps_g2_b_sqrt*mod_params->vdw_depth_cb_sqrt;
			vdw_erg = vdw_fn(b->g2, a->cb, rg2_b + mod_params->rcb,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 91.80) /* cg -- c */
		//	return erg;
        
		if (rg_a != 0.0 && mod_params->rc > 0.0) {
			depth = eps_g_a_sqrt*mod_params->vdw_depth_c_sqrt;
			vdw_erg = vdw_fn(a->g, b->c, rg_a + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rc > 0.0) {
			depth = eps_g2_a_sqrt*mod_params->vdw_depth_c_sqrt;
			vdw_erg = vdw_fn(a->g2, b->c, rg2_a + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rc > 0.0) {
			depth = eps_g_b_sqrt*mod_params->vdw_depth_c_sqrt;
			vdw_erg = vdw_fn(b->g, a->c, rg_b + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rc > 0.0) {
			depth = eps_g2_b_sqrt*mod_params->vdw_depth_c_sqrt;
			vdw_erg = vdw_fn(b->g2, a->c, rg2_b + mod_params->rc,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 88.06) /* cg -- n */
		//	return erg;
        
		if (rg_a != 0.0 && mod_params->rn > 0.0) {
			depth = eps_g_a_sqrt*mod_params->vdw_depth_n_sqrt;
			vdw_erg = vdw_fn(a->g, b->n, rg_a + mod_params->rn,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rn > 0.0) {
			depth = eps_g2_a_sqrt*mod_params->vdw_depth_n_sqrt;
			vdw_erg = vdw_fn(a->g2, b->n, rg2_a + mod_params->rn,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rn > 0.0) {
			depth = eps_g_b_sqrt*mod_params->vdw_depth_n_sqrt;
			vdw_erg = vdw_fn(b->g, a->n, rg_b + mod_params->rn,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rn > 0.0) {
			depth = eps_g2_b_sqrt*mod_params->vdw_depth_n_sqrt;
			vdw_erg = vdw_fn(b->g2, a->n, rg2_b + mod_params->rn,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
        
		//if (d2 > 67.33) /* cg - ca */
		//	return erg;
        
		if (rg_a != 0.0 && mod_params->rca > 0.0) {
			depth = eps_g_a_sqrt*mod_params->vdw_depth_ca_sqrt;
			vdw_erg = vdw_fn(a->g, b->ca, rg_a + mod_params->rca,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_a != 0.0 && mod_params->rca > 0.0) {
			depth = eps_g2_a_sqrt*mod_params->vdw_depth_ca_sqrt;
			vdw_erg = vdw_fn(a->g2, b->ca, rg2_a + mod_params->rca,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg_b != 0.0 && mod_params->rca > 0.0) {
			depth = eps_g_b_sqrt*mod_params->vdw_depth_ca_sqrt;
			vdw_erg = vdw_fn(b->g, a->ca, rg_b + mod_params->rca,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}
		if (rg2_b != 0.0 && mod_params->rca > 0.0) {
			depth = eps_g2_b_sqrt*mod_params->vdw_depth_ca_sqrt;
			vdw_erg = vdw_fn(b->g2, a->ca, rg2_b + mod_params->rca,depth, mod_params->rel_vdw_cutoff, mod_params->vdw_shift*depth, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
		}

	}

	/* 3. check even closer: this is the nonG-nonG case (CB-extended backbone) */
	if (d2 > backbone_constants[0]) /* o-o, farthest */
	/* 3.A  O-else */
		return erg;
    /* o - o, o - n, o - c */

	if (mod_params->ro > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->o, b->o, mod_params->ro + mod_params->ro,mod_params->vdw_depth_o_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_o_o, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	// only hard cutoff for N-O, the H-bond interactions will pull them in
	if (mod_params->ro > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->o, b->n, mod_params->ro + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->rn > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->n, b->o, mod_params->rn + mod_params->ro,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->ro > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->o, b->c, mod_params->ro + mod_params->rc,mod_params->vdw_depth_c_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_o, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->c, b->o, mod_params->rc + mod_params->ro,mod_params->vdw_depth_c_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_o, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}

   if ( d2 > backbone_constants[1])
		return erg; 

	/* 3.B  any-any, within N-N, mostly N with others */
	/* n - n, n - c, c - c, o - cb, o - ca, n - cb, c - cb */
	if (mod_params->rn > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->n, b->n, mod_params->rn + mod_params->rn,mod_params->vdw_depth_n_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_n_n, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	// only hard cutoff for N-C, the H-bond interactions will pull them in
	if (mod_params->rn > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->n, b->c, mod_params->rn + mod_params->rc,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->c, b->n, mod_params->rc + mod_params->rn,0, mod_params->rel_vdw_cutoff, 0, mod_params->vdw_clash_energy_at_hard_cutoff);
		erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->c, b->c, mod_params->rc + mod_params->rc,mod_params->vdw_depth_c_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_c, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->ro > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->o, b->cb, mod_params->ro + mod_params->rcb,mod_params->vdw_depth_cb_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_o, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->ro > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->o, mod_params->rcb + mod_params->ro,mod_params->vdw_depth_cb_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_o, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->ro > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->o, b->ca, mod_params->ro + mod_params->rca,mod_params->vdw_depth_ca_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_o, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->ro > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->o, mod_params->rca + mod_params->ro,mod_params->vdw_depth_ca_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_o, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rn > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->n, b->cb, mod_params->rn + mod_params->rcb,mod_params->vdw_depth_cb_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_cb, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rn > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->n, mod_params->rcb + mod_params->rn,mod_params->vdw_depth_cb_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_n, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->c, b->cb, mod_params->rc + mod_params->rcb,mod_params->vdw_depth_cb_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_c, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rc > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->c, mod_params->rcb + mod_params->rc,mod_params->vdw_depth_ca_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_ca, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	
	if ( d2 > backbone_constants[2])
		return erg; 
	
	/* 3.B  rest, within N-CA, mostly ca with all others */
	/* n - ca, c - ca, cb - ca, cb - cb, ca - ca */
	if (mod_params->rca > 0.0 && mod_params->rn > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->n, mod_params->rca + mod_params->rn,mod_params->vdw_depth_ca_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_n, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rn > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->n, b->ca, mod_params->rn + mod_params->rca,mod_params->vdw_depth_ca_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_n, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rc > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->c, b->ca, mod_params->rc + mod_params->rca,mod_params->vdw_depth_ca_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_c, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->rc > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->c, mod_params->rca + mod_params->rc,mod_params->vdw_depth_ca_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_c, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rcb > 0.0 && a->id != 'G' && b->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->cb, mod_params->rcb + mod_params->rcb,mod_params->vdw_depth_cb_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_cb, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->rcb > 0.0 && b->id != 'G') {
		vdw_erg = vdw_fn(a->ca, b->cb, mod_params->rca + mod_params->rcb,mod_params->vdw_depth_ca_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_cb, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rcb > 0.0 && mod_params->rca > 0.0 && a->id != 'G') {
		vdw_erg = vdw_fn(a->cb, b->ca, mod_params->rcb + mod_params->rca,mod_params->vdw_depth_ca_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_cb, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
	if (mod_params->rca > 0.0 && mod_params->rca > 0.0) {
		vdw_erg = vdw_fn(a->ca, b->ca, mod_params->rca + mod_params->rca,mod_params->vdw_depth_ca_ca, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_ca, mod_params->vdw_clash_energy_at_hard_cutoff);
			  erg += vdw_erg;
	}
    
	return erg;
}



/***********************************************************/
/****         CA-CA DISTANCE CUTOFF CALCULATION         ****/
/****               FOR  VDW INTERACTIONS               ****/
/***********************************************************/


double vdw_low_level(vector CA_1, vector CA_2, vector A_1, vector A_2, double radii_1, double radii_2){
  const double vdw_cutoff = 2.0;
  double y1 = sqrt(distance(CA_1,A_1));
  double y2 = sqrt(distance(CA_2,A_2));
  double ans = (y1 + y2 + vdw_cutoff*(radii_1 + radii_2));
  return ans*ans;	
}

/* This function returns the furthest CA_ - CA_ distance for backbone atom contacts.
   If the CA_ - CA_ distance of an amino acid pair is greater than the returned value,
   all backbone atom vdW contributions are 0.
   This function expects a peptide with first residue not glycine. */
double vdw_backbone_constants(Chain* chain, model_params *mod_params, FILE *outfile, int verbose){

  double bb_max = 0;
  double cut;

  if (chain->aa[1].id == 'G') stop("The vdW cutoff functions can only be calculated using a polypeptide chain not starting with glycine.");

  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].o,mod_params->ro,mod_params->ro);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"O-O %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].n,mod_params->ro,mod_params->rn);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"O-N %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].c,mod_params->ro,mod_params->rc);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"O-C %f\n",cut);
  //if (verbose) fprintf(stderr,"*******\n");
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].n,chain->aa[1].n,mod_params->rn,mod_params->rn);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"N-N %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].n,chain->aa[1].c,mod_params->rn,mod_params->rc); 
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"N-C %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].c,chain->aa[1].c,mod_params->rc,mod_params->rc); 
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"C-C %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].cb,mod_params->ro,mod_params->rcb);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"O-CB %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].ca,mod_params->ro,mod_params->rca);    
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"O-CA %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].n,chain->aa[1].cb,mod_params->rn,mod_params->rcb);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"N-CB %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].c,chain->aa[1].cb,mod_params->rc,mod_params->rcb);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"C-CB %f\n",cut);
  //if (verbose) fprintf(stderr,"*******\n");
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].n,chain->aa[1].ca,mod_params->rn,mod_params->rca); 
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"N-CA %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].c,chain->aa[1].ca,mod_params->rc,mod_params->rca); 
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"C-CA %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].cb,chain->aa[1].cb,mod_params->rcb,mod_params->rcb);
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"CB-CB %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].cb,chain->aa[1].ca,mod_params->rcb,mod_params->rca); 
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"CB-CA %f\n",cut);
  cut = vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].ca,mod_params->rca,mod_params->rca);	
  if (cut > bb_max) bb_max = cut;
  if (verbose) fprintf(stderr,"CA-CA %f\n",cut);
  //if (verbose) fprintf(stderr,"*******\n");

//  fprintf(outfile,"const double backbone_constants[3] = { %f, %f, %f };\n",vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].o,mod_params->ro,mod_params->ro),vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].n,mod_params->ro,mod_params->rn),vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].ca,mod_params->ro,mod_params->rca));
//  return vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].o,mod_params->ro,mod_params->ro);
  /* the largest of the 3 groups */
  if (verbose) fprintf(outfile,"const double backbone_constants[3] = { %f, %f, %f };\n",
		vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].o,mod_params->ro,mod_params->ro),
		vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].n,chain->aa[1].n,mod_params->rn,mod_params->rn),
		vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].n,chain->aa[1].ca,mod_params->rn,mod_params->rca));
  if (verbose) fprintf(outfile,"const double backbone_constants[3] = { %f, %f, %f };\n", bb_max, bb_max, bb_max);
  /* return the largest */
//  return vdw_low_level(chain->aa[1].ca,chain->aa[1].ca,chain->aa[1].o,chain->aa[1].o,mod_params->ro,mod_params->ro);

  mod_params->vdw_backbone_cutoff = bb_max;

  return bb_max;
}


/* This function returns the furthest CA_ - CA_ distance for [G__|G2_] - [G__|G2_] atom contacts.
   If the CA_ - CA_ distance of an amino acid pair is greater than the returned value,
   all [G__|G2_] - [G__|G2_] backbone atom vdW contributions are 0. */
double vdw_gamma_gamma(AA *a, AA *b, model_params *mod_params){

  double ans = mod_params->vdw_backbone_cutoff;

  if( a->etc & G2_ && b->etc & G2_){
	double temp = vdw_low_level(a->ca,b->ca,a->g2,b->g2,sidechain_vdw_radius(a->id,2, mod_params->sidechain_properties),sidechain_vdw_radius(b->id,2, mod_params->sidechain_properties));  
	//fprintf(stderr,"%c %d G2_ - %c %d G2_ is %f.\n",a->id,a->num,b->id,b->num,temp);
    if(temp > ans) {
	ans = temp;
    }
  }
  if( a->etc & G__ && b->etc & G2_){
	double temp = vdw_low_level(a->ca,b->ca,a->g,b->g2,sidechain_vdw_radius(a->id,1, mod_params->sidechain_properties),sidechain_vdw_radius(b->id,2, mod_params->sidechain_properties));  
    if(temp > ans) {
	ans = temp;
	//fprintf(stderr,"%c %d G__ - %c %d G2_ is %f.\n",a->id,a->num,b->id,b->num,temp);
    }
  }
  if( a->etc & G2_ && b->etc & G__){
	double temp = vdw_low_level(a->ca,b->ca,a->g2,b->g,sidechain_vdw_radius(a->id,2, mod_params->sidechain_properties),sidechain_vdw_radius(b->id,1, mod_params->sidechain_properties));  
    if(temp > ans) {
	ans = temp;
	//fprintf(stderr,"%c %d G2_ - %c %d G__ is %f.\n",a->id,a->num,b->id,b->num,temp);
    }
  }
  if( a->etc & G__ && b->etc & G__){
	double temp = vdw_low_level(a->ca,b->ca,a->g,b->g,sidechain_vdw_radius(a->id,1, mod_params->sidechain_properties),sidechain_vdw_radius(b->id,1, mod_params->sidechain_properties));  
    if(temp > ans) {
	ans = temp;
	//fprintf(stderr,"%c %d G__ - %c %d G__ is %f.\n",a->id,a->num,b->id,b->num,temp);
    }
  }
  return ans;	
} 

/* This function returns the furthest CA_ - CA_ distance for [G__|G2_] - ~[G__|G2_] atom contacts.
   If the CA_ - CA_ distance of an amino acid pair is greater than the returned value,
   all [G__|G2_] - ~[G__|G2_] vdW contributions are 0. */
double vdw_gamma_nongamma( AA *a,  AA *b, model_params *mod_params){

  double ans = mod_params->vdw_backbone_cutoff;

  if( a->etc & G2_){
	double temp = vdw_low_level(a->ca,b->ca,a->g2,b->o,sidechain_vdw_radius(a->id,2, mod_params->sidechain_properties),mod_params->ro);  
    if(temp > ans) ans = temp;
  }
  if( a->etc & G__ ){
	double temp = vdw_low_level(a->ca,b->ca,a->g,b->o,sidechain_vdw_radius(a->id,1, mod_params->sidechain_properties),mod_params->ro);  
    if(temp > ans) ans = temp;
  }
  if( b->etc & G2_){
	double temp = vdw_low_level(a->ca,b->ca,a->o,b->g2,mod_params->ro,sidechain_vdw_radius(b->id,2, mod_params->sidechain_properties));  
    if(temp > ans) ans = temp;
  }
  if(b->etc & G__){
	double temp = vdw_low_level(a->ca,b->ca,a->o,b->g,mod_params->ro,sidechain_vdw_radius(b->id,1, mod_params->sidechain_properties));  
    if(temp > ans) ans = temp;
  }
  return ans;	
}

/*This function expects a 26 residue peptide 
 * ABCDEFGHIGKLMNGPQRSTGVWGYZ and will output the
 * array vdw - this is not actually called in the program */	
void vdw_maxgamma_calc(Chain *chain, model_params *mod_params, FILE *outfile, int verbose) {

	int l, m;

	if (mod_params->use_gamma_atoms == NO_GAMMA) return;

	/* Reallocate memory to store the contact cutoff matrices */
	mod_params->vdw_gamma_gamma_cutoff = (double *)realloc(mod_params->vdw_gamma_gamma_cutoff, 702 * sizeof(double));
	mod_params->vdw_gamma_nongamma_cutoff = (double *)realloc(mod_params->vdw_gamma_nongamma_cutoff, 702 * sizeof(double));

	/* Quick check of side chains */
	for (int i=1; i<chain->NAA; i++) {
		if ( ((chain->aa[i].id == 'I') || (chain->aa[i].id == 'T') || (chain->aa[i].id == 'V')) &&
		     ((chain->aa[i].etc & G2_) == 0) ) stop("I, T and V must have an explicit G2_ atom for the vdW cutoff matrix generation.");
		if ( ((chain->aa[i].id != 'G') && (chain->aa[i].id != 'A') ) &&
		     ((chain->aa[i].etc & G__) == 0) ) stop("all residues that are not G or S must have an explicit G__ atom for the vdW cutoff matrix generation.");
	}
	if (verbose) {
		fprintf(stderr,"The polypeptide chain on which the vdW cutoff matrices are calculated:\n");
		fprintf(stderr,"aa   G__   G2_\n");
		fprintf(stderr,"==============\n");
		for (int i=1; i<chain->NAA; i++) {
			fprintf(stderr,"%c     %d     %d\n", chain->aa[i].id, ((chain->aa[i].etc & G__) != 0), ((chain->aa[i].etc & G2_) != 0));
		}
	}

	/* Calculate the vdW cutoff matrices */
	for(m = 1; m <= 26; m++) { /* blank row, not a valid amino acid */
		mod_params->vdw_gamma_gamma_cutoff[m-1] = mod_params->vdw_backbone_cutoff;
		mod_params->vdw_gamma_nongamma_cutoff[m-1] = mod_params->vdw_backbone_cutoff;
	}
	/* The 26x26 array of gamma - gamma contact cutoffs */
	for(l = 1; l <= 26; l++){
	    for(m = 1; m <= 26; m++){
		mod_params->vdw_gamma_gamma_cutoff[(l)*26+(m-1)] = vdw_gamma_gamma(&(chain->aa[l]),&(chain->aa[m]),mod_params);
	    }
	}
	/* The 26x26 array of gamma - nongamma contact cutoffs */
	for(l = 1; l <= 26; l++){
	  for(m = 1; m <= 26; m++){
		mod_params->vdw_gamma_nongamma_cutoff[(l)*26+(m-1)] = vdw_gamma_nongamma(&(chain->aa[l]),&(chain->aa[m]),mod_params);
	  }
	}

}

/* This function calculates the vdW cutoff matrices
   for the model parameters used.
 * ABCDEFGHIGKLMNGPQRSTGVWGYZ and will output the
 * array vdw - this is not actually called in the program */	
void vdw_cutoff_distances_calculate(simulation_params *sim_params, FILE *outfile, int verbose) {

	/* build test peptide */
	simulation_params *my_sim_params = (simulation_params *)malloc(sizeof(simulation_params));
	sim_params_copy(my_sim_params,sim_params);
	my_sim_params->infile = NULL;
	my_sim_params->outfile = NULL;
	my_sim_params->checkpoint_file = NULL;
	Chain *chain = (Chain *)malloc(sizeof(Chain)*28); chain->NAA = 0;
        Chaint* chaint = (Chaint *)malloc(sizeof(Chaint));
      	chain->NAA =0; chain->Nchains = 0;
      	chain->aa = NULL; chain->xaa = NULL; chain->xaa_prev = NULL; chain->erg = NULL;
      	chaint->aat = NULL; chaint->xaat = NULL; chaint->xaat_prev = NULL; chaint->ergt = NULL;  
	if (my_sim_params->seq) free(my_sim_params->seq);
	copy_string(&(my_sim_params->seq),"ABCDEFGHIGKLMNGPQRSTGVWGYZ");
	build_peptide_from_sequence(chain,chaint,my_sim_params->seq, my_sim_params);
	param_finalise(my_sim_params);
	free(my_sim_params);

	/* vdW cutoff for backbone atoms */
	/* Use extended cutoff if requested. */
	if (sim_params->protein_model.vdw_use_extended_cutoff) {
		sim_params->protein_model.vdw_backbone_cutoff = sim_params->protein_model.vdw_extended_cutoff;
	} else {
		vdw_backbone_constants(chain,&(sim_params->protein_model), stderr, /* verbose = */ 0);
	}

	/* gamma-gamma and gamma-nongamma vdW cutoffs */
	vdw_maxgamma_calc(chain, &(sim_params->protein_model), stderr, /* verbose = */ 0);
	/* to be on the safe side, let's add 1 to all */
	for (int i=0; i<702 ; i++) {
		sim_params->protein_model.vdw_gamma_gamma_cutoff[i] += 10.0;
		sim_params->protein_model.vdw_gamma_nongamma_cutoff[i] += 10.0;
	}

	if (verbose) print_vdw_cutoff_distances(&(sim_params->protein_model),outfile);

	/* free memory */
	freemem_chaint(chaint);
	free(chaint);
	freemem_chain(chain); //free amino acid chain and energy matrix
	free(chain);

}

/* Updating NAA (number of amino acids),
            seq (amino acid sequence),
            sequence (amino acid sequence with '_' separators for chain breaks) and
            N (number of NS active points)
   of sim params from chain */
void update_sim_params_from_chain(Chain *chain,simulation_params *sim_params) {

  //TODO: update sim_params->N as well
  sim_params->N = 0;
  //fprintf(stderr,"updating sim_param NAA to %d and then seq.\n",chain->NAA);
  sim_params->NAA = chain->NAA;
  sim_params->Nchains = chain->Nchains;
  sim_params->seq = realloc(sim_params->seq,(sim_params->NAA+1) * sizeof(char));
  // The first letter is A for some weird reason, it does not code any amino acid.
  sim_params->seq[0] = 'A';
  for (int i=1; i<sim_params->NAA; i++) {
    sim_params->seq[i] = chain->aa[i].id;
  }
  sim_params->seq[sim_params->NAA] = '\0';

  /* sequence for multi-chain proteins */
  sim_params->sequence = realloc(sim_params->sequence,(sim_params->NAA+sim_params->Nchains) * sizeof(char));
  // The first letter is A for some weird reason, it does not code any amino acid.
//  fprintf(stderr,"A%c",chain->aa[1].id);
  sim_params->sequence[0] = 'A';
  sim_params->sequence[1] = chain->aa[1].id;
  int next = 2;
  int Nchains = 1;
  for (int i=2; i<sim_params->NAA; i++) {
//    fprintf(stderr,"%c",chain->aa[i].id);
    /* add chain break if necessary */
    if (chain->aa[i].chainid != chain->aa[i-1].chainid){
//      fprintf(stderr,"gap found between i = %d (%d) and i-1 = %d (%d)\n",i,chain->aa[i].chainid,i-1,chain->aa[i-1].chainid);
      sim_params->sequence[next++] = '_';
      Nchains ++;
    }
    sim_params->sequence[next++] = chain->aa[i].id;
  }
  sim_params->sequence[next] = '\0';
//  fprintf(stderr,"\n");
//  fprintf(stderr,"\nUpdated sim_params->sequence: %s (%d)\n",sim_params->sequence,next);
  if (Nchains != chain->aa[sim_params->NAA-1].chainid) {
	fprintf(stderr,"Nchains = %d, last_chainid = %d\n",Nchains,chain->aa[sim_params->NAA-1].chainid);
  	stop("vdw.c: The number of chains != the last chain ID.\n");
  }
  if (Nchains != chain->Nchains) {
	fprintf(stderr,"Nchains = %d, nchains = %d\n",Nchains,chain->Nchains);
	stop("vdw.c: The number of chains != the last chain ID.\n");
  }
  sim_params->Nchains = Nchains;
}
