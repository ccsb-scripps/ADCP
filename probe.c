/*
** These routines calculate and print important polypeptide characteristics.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

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
#include"flex.h"

#ifdef PARALLEL
#include<mpi.h>
#endif
#include"checkpoint_io.h"

#define Distb(I,J)     biasmap->distb[(I) * biasmap->NAA + (J)]

/* test registry */
struct TEST {
	void (*func) (Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
	char *name;
	int init_mask;
};

/* number of tests should not exceed the integer bits (32) */
struct TEST test[] = {
	{rgyr, "Centroid and radius of gyration (default)", 0x01 }, /* before init */
	{pdbout, "Snapshots in PDB format (default)", 0x11 }, /* ((before and)) after init */
	{phipsi, "Ramachandran dihedral angles and backbone angle and side chain dihedral angles", 0x10 }, /* before and after init */
	{sstructure, "Secondary structure assignment from dihedral angles", 0x10 }, /* after init */
	{n_contacts, "Total number of contacts", 0x1 }, /* before init */
	{cm_txt, "Text contact map !FOR SIMULATION INPUT GENERATION USE cm_ideal INSTEAD!", 0x1 }, /* before init */
	{cm_pbm, "Contact map image in PBM format", 0x1 }, /* before init */
	{hbtot, "Total number of hydrogen bonds", 0x1 }, /* before init */
	{hpattern, "Local binary H-bonding donor and acceptor patterns", 0x1 }, /* before init */
	{hbss, "DSSP-like secondary structure assignment", 0x1 }, /* before init */
	{writhe, "Writhe of alpha-carbon trace", 0x1 }, /* before init */
	{ergtot, "Total and local energy (default only for NS)", 0x10 }, /* after init */
	{stepsize, "Crankshaft amplitude, acceptance rate, other diagnostics", 0x10 }, /* after init */
	{energy_gradient_wrt_parameters, "The gradient of the energy with respect to its parameters, calculated by finite differences", 0x10 }, /* only after init */
	{n_native_contacts, "Number of native contacts", 0x10 },
	{evidence,"log(Evidence) estimate (only NS, default)", 0x11 },
	{information,"Information estimate (only NS, default)", 0x11 },
	{cm_ideal, "Ideal contact map", 0x11 },
	//{cm_ideal_4, "Ideal contact map including i,i+4 alpha contacts", 0x11 },
	//{cm_native_go, "Native contact map including all Cbeta-Cbeta contacts (useful for the calculation of the number of native contacts).", 0x01 },
	{test_flex,"Test Flex",0x10 }, /* only after init */
	{energy_contributions, "Energy contributions of all different terms", 0x10 }, /* only after init */
	{exclude_energy_contributions, "Exclude energy of amino acid pairs", 0x10 }, /* only after init */
	{cos_dihedral_naac, "cos of N(i)-CA(i)-CA(i+1)-C(i+1) dihedral angle", 0x11 }, /* before and after init */
	{strand_bias_distances, "CB-CB and CA-CA distances of beta contacts", 0x11 }, /* before and after init */
	{hydrophobic_distances, "CB-CB and CA-CA distances of HH and HP contacts", 0x10 }, /* after init */
	{number_of_contacts, "Number of contacts", 0x11 },
	{initialize_displacement, "The displacement of atoms and generated clashes during initialization", 0x1 }, /* only before init */
	{hydrophobic_rgyr, "Centroid and radius of gyration of only the hydrophobic residues", 0x01 }, /* before init */
	{CA_geometry, "Various geometrical parameters using CA atoms (atomic distances, angles, dihedrals)", 0x01 }, /* before init */
	{atomic_distances,"List all atomic distances, useful for LJ testing",0x10 },
	{hbond_geometry, "All geometrical parameters of Hbonds: d(O,H), <(C,O,H), <(O,H,N)",0x10 },
	{vdw_max_gamma, "Calculate vdW max gamma-(non)gamma tables for the given vdW parameters. Use peptide ABCDEFGHIGKLMNGPQRSTGVWGYZ, and do not forget to specify the VDW parameters.",0x10 },
//	{vdw_contributions, "Calculate vdW contributions between all atoms.",0x10 },
	//{hbond_pattern, "Calculate H-bond pattern.",0x01 },
	{checkpoint_out, "Snapshots in CHK format", 0x11 }, /* ((before and)) after init */
//	{cm_alpha_8, "Contact map with contacts where d(CA-CA) < 8 Angstroms.", 0x01 },
	//{fasta, "Print fasta sequence.", 0x01 },
	{NULL, NULL, 0x0 }
};

/* helper buffers */
//char *ss = NULL;
int *hbm = NULL;
//char *map = NULL;
#define Hbm(I,J)     hbm[(I) * chain->NAA + (J)]
#define Map(I,J)     map[(I) * chain->NAA + (J)]

void evidence(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm){
; //see nested.c	
}

void information(Chain *chain,Biasmap *biasmap,  simulation_params *sim_params, void *mpi_comm){
; //see nested.c	
}

void tests(Chain *chain, Biasmap *biasmap, unsigned int t, simulation_params *sim_params, int init_mask, void *mpi_comm)
{
	int i;

	for (i = 0; test[i].func != NULL && t; i++, t >>= 1)
		//do test if its mask is given and if it is meant to be calculated (some are only valid before/after initialization)
		if ((t & 0x1) && (init_mask & test[i].init_mask))
			(test[i].func) (chain,biasmap,sim_params, mpi_comm);
}

void helps(void)
{
	int i;
	unsigned int t = 0x1;

	fprintf(stderr, "TestMASK Index:\n");

	for (i = 0; test[i].func != NULL && t; i++, t <<= 1)
		fprintf(stderr, "%8x %s\n", t, test[i].name);
}

void test_flex(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm){
  if(sim_params->flex_params.number_of_processors > 0){
    Chain *input_chains;
    double *rmsd;
    initialize_flex(chain,&input_chains,biasmap,sim_params,&rmsd);
    output_and_run_flex(chain,biasmap,input_chains,sim_params,rmsd);
    finalize_flex(chain,&input_chains,sim_params,&rmsd);
  }
}


/* various geometrical parameters using CA atoms*/
void CA_geometry(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm) {

//   for (i=0; i<sim_params->geometry.nbonds; i++) {
//      int aa1 = sim_params->geometry.bond1[i];
//      int aa2 = sim_params->geometry.bond2[i];
//      fprintf(sim_params->outfile,"distance %d %d %g", sqrt(distance(chain->aa[aa1].ca,chain->aa[aa2].ca)));
//   }
      fprintf(sim_params->outfile,"distance  1 16 %g", sqrt(distance(chain->aa[1].ca,chain->aa[16].ca)));
      fprintf(sim_params->outfile,"distance  3 14 %g", sqrt(distance(chain->aa[3].ca,chain->aa[14].ca)));
      fprintf(sim_params->outfile,"distance  5 12 %g", sqrt(distance(chain->aa[5].ca,chain->aa[12].ca)));
      fprintf(sim_params->outfile,"distance  7 10 %g", sqrt(distance(chain->aa[7].ca,chain->aa[10].ca)));

      if (chain->NAA>17) {
        fprintf(sim_params->outfile,"distance  18 15 %g", sqrt(distance(chain->aa[18].ca,chain->aa[15].ca)));
        fprintf(sim_params->outfile,"distance  20 13 %g", sqrt(distance(chain->aa[20].ca,chain->aa[13].ca)));
        fprintf(sim_params->outfile,"distance  22 11 %g", sqrt(distance(chain->aa[22].ca,chain->aa[11].ca)));
        fprintf(sim_params->outfile,"distance  24  9 %g", sqrt(distance(chain->aa[24].ca,chain->aa[9].ca)));
      }
}

/* various geometrical parameters using CA atoms*/
void vdw_contributions(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm) {

	model_params *mod_params = &sim_params->protein_model;

	//Assign function pointer to the vdW model
	double (*vdw_fn) (vector r1, vector r2, double Rmin, double depth, double vdw_rel_cutoff, double energy_shift, double clash_at_vdw_cutoff) = NULL;
	if (mod_params->vdw_potential == HARD_CUTOFF_VDW_POTENTIAL) {
		vdw_fn = vdw_hard_cutoff;
	} else if (mod_params->vdw_potential == LJ_VDW_POTENTIAL) {
		vdw_fn = vdw_lj;
	} else {
		stop("Clash cannot be calculated without a valid vdW potential.");
	}


	for (int i=1; i<chain->NAA; i++) {
	    for (int j=i; j<chain->NAA; j++) {
		fprintf(sim_params->outfile,"%d %d interactions:",i,j);
		AA *a = &(chain->aa[i]);
		AA *b = &(chain->aa[j]);
		fprintf(sim_params->outfile,"CA_ CA_ %g\n",vdw_fn(a->ca, b->ca, mod_params->rca + mod_params->rca,mod_params->vdw_depth_ca_ca, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_ca, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CA_ CB_ %g\n",vdw_fn(a->ca, b->cb, mod_params->rca + mod_params->rcb,mod_params->vdw_depth_ca_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_cb, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CA_ C__ %g\n",vdw_fn(a->ca, b->c, mod_params->rca + mod_params->rc,mod_params->vdw_depth_ca_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_c, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CA_ N__ %g\n",vdw_fn(a->ca, b->n, mod_params->rca + mod_params->rn,mod_params->vdw_depth_ca_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_n, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CA_ O__ %g\n",vdw_fn(a->ca, b->o, mod_params->rca + mod_params->ro,mod_params->vdw_depth_ca_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CB_ CA_ %g\n",vdw_fn(a->cb, b->ca, mod_params->rca + mod_params->rca,mod_params->vdw_depth_ca_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_cb, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CB_ CB_ %g\n",vdw_fn(a->cb, b->cb, mod_params->rcb + mod_params->rcb,mod_params->vdw_depth_cb_cb, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_cb, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CB_ C__ %g\n",vdw_fn(a->cb, b->c, mod_params->rcb + mod_params->rc,mod_params->vdw_depth_cb_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_c, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CB_ N__ %g\n",vdw_fn(a->cb, b->n, mod_params->rcb + mod_params->rn,mod_params->vdw_depth_cb_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_n, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"CB_ O__ %g\n",vdw_fn(a->cb, b->o, mod_params->rcb + mod_params->ro,mod_params->vdw_depth_cb_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"C__ CA_ %g\n",vdw_fn(a->c, b->ca, mod_params->rcb + mod_params->rca,mod_params->vdw_depth_ca_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_c, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"C__ CB_ %g\n",vdw_fn(a->c, b->cb, mod_params->rc + mod_params->rcb,mod_params->vdw_depth_cb_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_c, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"C__ C__ %g\n",vdw_fn(a->c, b->c, mod_params->rc + mod_params->rc,mod_params->vdw_depth_c_c, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_c, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"C__ N__ %g\n",vdw_fn(a->c, b->n, mod_params->rc + mod_params->rn,mod_params->vdw_depth_c_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_n, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"C__ O__ %g\n",vdw_fn(a->c, b->o, mod_params->rc + mod_params->ro,mod_params->vdw_depth_c_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"N__ CA_ %g\n",vdw_fn(a->n, b->ca, mod_params->rn + mod_params->rca,mod_params->vdw_depth_ca_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_n, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"N__ CB_ %g\n",vdw_fn(a->n, b->cb, mod_params->rn + mod_params->rcb,mod_params->vdw_depth_cb_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_n, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"N__ C__ %g\n",vdw_fn(a->n, b->c, mod_params->rn + mod_params->rc,mod_params->vdw_depth_c_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_n, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"N__ N__ %g\n",vdw_fn(a->n, b->n, mod_params->rn + mod_params->rn,mod_params->vdw_depth_n_n, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_n_n, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"N__ O__ %g\n",vdw_fn(a->n, b->o, mod_params->rn + mod_params->ro,mod_params->vdw_depth_n_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_n_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"O__ CA_ %g\n",vdw_fn(a->o, b->ca, mod_params->ro + mod_params->rca,mod_params->vdw_depth_ca_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_ca_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"O__ CB_ %g\n",vdw_fn(a->o, b->cb, mod_params->ro + mod_params->rcb,mod_params->vdw_depth_cb_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_cb_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"O__ C__ %g\n",vdw_fn(a->o, b->c, mod_params->ro + mod_params->rc,mod_params->vdw_depth_c_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_c_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"O__ N__ %g\n",vdw_fn(a->o, b->n, mod_params->ro + mod_params->rn,mod_params->vdw_depth_n_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_n_o, mod_params->vdw_clash_energy_at_hard_cutoff));
		fprintf(sim_params->outfile,"O__ O__ %g\n",vdw_fn(a->o, b->o, mod_params->ro + mod_params->ro,mod_params->vdw_depth_o_o, mod_params->rel_vdw_cutoff, mod_params->vdw_Eshift_o_o, mod_params->vdw_clash_energy_at_hard_cutoff));
	    }
	}
}

/* the displacement of atoms and generated clashes during initialization */
void initialize_displacement(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm) {

	Chain chain_init;
	chain_init.NAA = 0;
	chain_init.aa = NULL; chain_init.xaa = NULL; chain_init.erg = NULL; chain_init.xaa_prev = NULL;
        allocmem_chain(&chain_init,chain->NAA,chain->Nchains);
	Chaint chaint;
	chaint.aat = NULL; chaint.xaat = NULL; chaint.ergt = NULL; chaint.xaat_prev = NULL;
	Chain displacements;
	displacements.NAA = 0;
	displacements.aa = NULL; displacements.xaa = NULL; displacements.erg = NULL; displacements.xaa_prev = NULL;
        allocmem_chain(&displacements,chain->NAA,chain->Nchains);
	
	/* calculate initialized distances */
	copybetween(&chain_init,chain);
//	fprintf(sim_params->outfile,"before init chain_init\n");
	chkpeptide(chain_init.aa, chain_init.NAA, &(sim_params->protein_model));
//	pdbprint(chain_init.aa, chain_init.NAA, &(sim_params->protein_model), sim_params->outfile, NULL);
////	energy_init(&chain_init,biasmap,&(sim_params->protein_model));
	initialize(&chain_init,&chaint,sim_params); // PEPTIDE MODIFICATION!!
////	energy_matrix_init(chain,biasmap,&(sim_params->protein_model));
//	fprintf(sim_params->outfile,"after init chain_init\n");
//	pdbprint(chain_init.aa, chain_init.NAA, &(sim_params->protein_model), sim_params->outfile, NULL);

	double threshold = 0.5;
	double dist;
	fprintf(stderr,"Listing displacements larger than %g Angstroms (not listing G__ atoms):\n",threshold);
	/* calculate displacements of all atoms */
	for (int i=1; i<chain_init.NAA; i++) {
		if (chain_init.aa[i].etc & CA_) { /* CA_ */
			dist = sqrt(distance(chain_init.aa[i].ca, chain->aa[i].ca));
			displacements.aa[i].ca[0] = dist;
			if (dist > threshold) fprintf(stderr,"Big displacement of %c%d CA_: %g\n",chain_init.aa[i].id,chain_init.aa[i].num,dist);
		}
		if (chain_init.aa[i].etc & C__) { /* C__ */
			dist = sqrt(distance(chain_init.aa[i].c, chain->aa[i].c));
			displacements.aa[i].c[0] = dist;
			if (dist > threshold) fprintf(stderr,"Big displacement of %c%d C__: %g\n",chain_init.aa[i].id,chain_init.aa[i].num,dist);
		}
		if (chain_init.aa[i].etc & N__) { /* N_ */
			dist = sqrt(distance(chain_init.aa[i].n, chain->aa[i].n));
			displacements.aa[i].n[0] = dist;
			if (dist > threshold) fprintf(stderr,"Big displacement of %c%d N__: %g\n",chain_init.aa[i].id,chain_init.aa[i].num,dist);
		}
		if (chain_init.aa[i].etc & O__) { /* O__ */
			dist = sqrt(distance(chain_init.aa[i].o, chain->aa[i].o));
			displacements.aa[i].o[0] = dist;
			if (dist > threshold) fprintf(stderr,"Big displacement of %c%d O__: %g\n",chain_init.aa[i].id,chain_init.aa[i].num,dist);
		}
		if (chain_init.aa[i].etc & CB_) { /* CB_ */
			dist = sqrt(distance(chain_init.aa[i].cb, chain->aa[i].cb));
			displacements.aa[i].cb[0] = dist;
			if (dist > threshold) {
				fprintf(stderr,"Big displacement of %c%d CB_: %g\n",chain_init.aa[i].id,chain_init.aa[i].num,dist);
				//fprintf(stderr," old %g %g %g and new %g %g %g\n",chain->aa[i].cb[0],chain->aa[i].cb[1],chain->aa[i].cb[2],chain_init.aa[i].cb[0],chain_init.aa[i].cb[1],chain_init.aa[i].cb[2]);
			}
		}
		/* gamma atoms will always move a lot -- don't list them */
		if (chain_init.aa[i].etc & G__) { /* G__ */
			dist = sqrt(distance(chain_init.aa[i].g, chain->aa[i].g));
			displacements.aa[i].g[0] = dist;
		//	if (dist > threshold) fprintf(stderr,"Big displacement of %c%d G__: %g\n",chain_init.aa[i].id,chain_init.aa[i].num,dist);
		}
		if (chain_init.aa[i].etc & G2_) { /* G2_ */
			dist = sqrt(distance(chain_init.aa[i].g2, chain->aa[i].g2));
			displacements.aa[i].g2[0] = dist;
		//	if (dist > threshold) fprintf(stderr,"Big displacement of %c%d G2_: %g\n",chain_init.aa[i].id,chain_init.aa[i].num,dist);
		}
	}

	int *G_to_delete = malloc(chain_init.NAA * sizeof(int));
	int *G2_to_delete = malloc(chain_init.NAA * sizeof(int));
	for (int i=0; i<chain_init.NAA; i++) {
	   G_to_delete[i] = 0;
	   G2_to_delete[i] = 0;
	   //fprintf(sim_params->outfile,"%x ",chain_init.aa[i].etc);
	}
	//fprintf(sim_params->outfile,"\n");
	
	/* calculate exclude energies with and without each gamma atom */
	for (int i=1; i<chain_init.NAA; i++) {
	   /* G__ */
	   if (chain_init.aa[i].etc & G__) {
		/* exclude */
		double exclude_with = 0;
		for (int j=1; j<chain_init.NAA; j++) {
		    if (abs(i-j)<2) continue;
		    exclude_with += exclude(&(chain_init.aa[i]), &(chain_init.aa[j]), 0.0, &(sim_params->protein_model));
		}
		double exclude_without = 0;
		chain_init.aa[i].etc &= ~G__;
		for (int j=1; j<chain_init.NAA; j++) {
		    if (abs(i-j)<2) continue;
		    exclude_without += exclude(&(chain_init.aa[i]), &(chain_init.aa[j]), 0.0, &(sim_params->protein_model));
		}
		chain_init.aa[i].etc |= G__;
		if (exclude_with > exclude_without + 5.0) {
//			fprintf(stderr,"XXX ");
			G_to_delete[i] = 1;
		}
//		fprintf(stderr,"The clash energy of %c%d with and without G__ atom is %g and %g\n",
//			chain_init.aa[i].id, chain_init.aa[i].num,
//			exclude_with, exclude_without);
		/* clash */
		double clash_with = clash(&(chain_init.aa[i]), &(sim_params->protein_model));
		chain_init.aa[i].etc &= ~G__;
		double clash_without = clash(&(chain_init.aa[i]), &(sim_params->protein_model));
		chain_init.aa[i].etc |= G__;
		if (clash_with > clash_without + 5.0) {
			G_to_delete[i] = 1;
		fprintf(stderr,"The clash energy of %c%d with and without G__ atom is %g and %g\n",
			chain_init.aa[i].id, chain_init.aa[i].num,
			clash_with, clash_without);
		}
	   }
	   /* G2_ */
	   if (chain_init.aa[i].etc & G2_) {
		/* exclude */
		double exclude_with = 0;
		for (int j=1; j<chain_init.NAA; j++) {
		    if (abs(i-j)<2) continue;
		    exclude_with += exclude(&(chain_init.aa[i]), &(chain_init.aa[j]), 0.0, &(sim_params->protein_model));
		}
		double exclude_without = 0;
		chain_init.aa[i].etc &= ~G2_;
		for (int j=1; j<chain_init.NAA; j++) {
		    if (abs(i-j)<2) continue;
		    exclude_without += exclude(&(chain_init.aa[i]), &(chain_init.aa[j]), 0.0, &(sim_params->protein_model));
		}
		chain_init.aa[i].etc |= G2_;
		if (exclude_with > exclude_without + 5.0) {
//			fprintf(stderr,"XXX ");
			G2_to_delete[i] = 1;
		}
//		fprintf(stderr,"The clash energy of %c%d with and without G2_ atom is %g and %g\n",
//			chain_init.aa[i].id, chain_init.aa[i].num,
//			exclude_with, exclude_without);
		/* clash */
		double clash_with = clash(&(chain_init.aa[i]), &(sim_params->protein_model));
		chain_init.aa[i].etc &= ~G2_;
		double clash_without = clash(&(chain_init.aa[i]), &(sim_params->protein_model));
		chain_init.aa[i].etc |= G2_;
		if (clash_with > clash_without + 5.0) {
			G_to_delete[i] = 1;
		fprintf(stderr,"The clash energy of %c%d with and without G2_ atom is %g and %g\n",
			chain_init.aa[i].id, chain_init.aa[i].num,
			clash_with, clash_without);
		}
	   }
	}

	/* remove clashing gamma atoms */
	fprintf(sim_params->outfile,"Removed clashing gamma atoms: ");
	for (int i=1; i<chain_init.NAA; i++) {
	   if (G_to_delete[i]) {
		fprintf(sim_params->outfile,"%c%d-G ",chain_init.aa[i].id,i);
		chain_init.aa[i].etc &= ~G__;
	   }
	   if (G2_to_delete[i]) {
		fprintf(sim_params->outfile,"%c%d-G2 ",chain_init.aa[i].id,i);
		chain_init.aa[i].etc &= ~G2_;
	   }
	   //fprintf(sim_params->outfile,"%x ",chain_init.aa[i].etc);
	}
	fprintf(sim_params->outfile,"\n");

	/* check that there are no more clashes */
	for (int i=1; i<chain_init.NAA; i++) {
		double exclude_all = 0;
		for (int j=1; j<chain_init.NAA; j++) {
		    if (abs(i-j)<2) continue;
		    exclude_all += exclude(&(chain_init.aa[i]), &(chain_init.aa[j]), 0.0, &(sim_params->protein_model));
		}
		if (exclude_all > 5.0)
		    fprintf(stderr,"WARNING! The clash energy of %c%dafter removing clashing gamma atoms is %g\n",
			chain_init.aa[i].id, chain_init.aa[i].num, exclude_all);
	}

	/* print PDB */
	fprintf(sim_params->outfile,"--+ PDB without clashing gamma atoms +--\n");
	fprintf(stderr,"Outputting initialized PDB?\n");
	pdbprint(chain_init.aa, chain_init.NAA, &(sim_params->protein_model), sim_params->outfile, NULL);

//	/* calculate clashes and print displacements and vdW energy of clashing atoms */
//	for (int i=1; i<chain_init.NAA; i++) {
//	   for (int j=i+1; j<chain_init.NAA; j++) {
//		double vdw_energy;
//		/* CA_ - CA_ */
//		if (j!=i+1){
//		    vdw_energy = vdw(chain_init.aa[i].ca,chain_init.aa[j].ca,sim_params->protein_model.rca+sim_params->protein_model.rca,sim_params->protein_model.vdw_depth_ca_ca,sim_params->protein_model.rel_vdw_cutoff,sim_params->protein_model.vdw_Eshift_ca_ca)>0.0;
//		    if (vdw_energy > 0.0) {
//			fprintf(stderr,"Clash %c%d CA_ and %c%d CA_ with displacements %g and %g and distance %g and energy %g, \n",
//				chain_init.aa[i].id,chain_init.aa[i].num,
//				chain_init.aa[j].id,chain_init.aa[j].num,
//				displacements.aa[i].ca[0],displacements.aa[j].ca[0],
//				sqrt(distance(chain_init.aa[i].ca,chain_init.aa[j].ca)), vdw_energy);
//		    }
//		    /* all */
//		    vdw_energy += exclude(&(chain_init.aa[i]), &(chain_init.aa[j]), 0.0, &(sim_params->protein_model));
//		    if (vdw_energy > 0.0) {
//			fprintf(stderr,"Clash %c%d and %c%d with energy %g, \n",
//				chain_init.aa[i].id,chain_init.aa[i].num,
//				chain_init.aa[j].id,chain_init.aa[j].num,
//				vdw_energy);
//		    }
//		    if (vdw_energy > 5.0) {
//			fprintf(stderr,"Clash %c%d and %c%d with energy %g",
//				chain_init.aa[i].id,chain_init.aa[i].num,
//				chain_init.aa[j].id,chain_init.aa[j].num,
//				vdw_energy);
//			if (chain_init.aa[i].etc & G__) fprintf(stderr,", 1st G__ displacement %g",displacements.aa[i].g[0]);
//			if (chain_init.aa[i].etc & G2_) fprintf(stderr,", 1st G2_ displacement %g",displacements.aa[i].g2[0]);
//			if (chain_init.aa[j].etc & G__) fprintf(stderr,", 2nd G__ displacement %g",displacements.aa[j].g[0]);
//			if (chain_init.aa[j].etc & G2_) fprintf(stderr,", 2nd G2_ displacement %g\n",displacements.aa[j].g2[0]);
//		    }
//		}
//	   }
//	}
//

        freemem_chain(&chain_init);

}

void torsion(Chain *chain, Biasmap *biasmap,double *phi, double *psi, int i)
{

	vector nca, cac;	/* N-Ca and Ca-C bonds */

	subtract(nca, chain->aa[i].ca, chain->aa[i].n);
	subtract(cac, chain->aa[i].c, chain->aa[i].ca);

	if (chain->aa[i].chainid != chain->aa[i-1].chainid) {
		*phi = M_180_PI * dihedral(chain->xaa_prev[chain->aa[i].chainid][1], nca, cac);
	} else {
		*phi = M_180_PI * dihedral(chain->xaa[i - 1][1], nca, cac);
	}
	*psi = M_180_PI * dihedral(nca, cac, chain->xaa[i][1]);
}

void all_torsions(Chain *chain,Biasmap *biasmap, double *phi, double *psi, double *chi1, double *chi2, int i)
{
	vector nca, cacb, cbg, cbg2;	/* N-Ca, Ca-Cb, Cb-G, Cb-G2 bonds */
	vector cbca;
	double NaN;
	NaN = strtod("NaN", NULL);

	/* phi, psi */
	torsion(chain,biasmap,phi, psi, i);

	/* chi1, chi2 */
	*chi1 = NaN;
	*chi2 = NaN;
	if ( chain->aa[i].etc & G__ ) { /* at least one gamma atom */
	   subtract(nca, chain->aa[i].ca, chain->aa[i].n);
	   subtract(cacb, chain->aa[i].cb, chain->aa[i].ca);
	   subtract(cbca, chain->aa[i].ca, chain->aa[i].cb);
	   subtract(cbg, chain->aa[i].g, chain->aa[i].cb);
	   *chi1 = M_180_PI * dihedral(nca,cacb,cbg);
	   if ( chain->aa[i].etc & G2_ ) { /* two gamma atoms */
	      subtract(cbg2, chain->aa[i].g2, chain->aa[i].cb);
	      *chi2 = M_180_PI * dihedral(nca,cacb,cbg2);
	   }
	} else { /* first gamma atom is not present */
	   if ( chain->aa[i].etc & G2_ ) { /* only second gamma atoms */
	      fprintf(stderr,"WARNING: amino acid %c%d has the second gamma atom, but not the first one!",chain->aa[i].id,i);
	      subtract(cbg2, chain->aa[i].g2, chain->aa[i].cb);
	      *chi2 = M_180_PI * dihedral(nca,cacb,cbg2);
	   }
	}

}

void sstructure(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, ch;
	double phi, psi;

	for (i = 1; i < chain->NAA; i++) {
		torsion(chain,biasmap,&phi, &psi, i);

		if (!(0. < phi && phi < 120.) && (-120. < psi && psi < 60.))
			ch = 'H';
		else
			ch = 'C';
		fputc(ch,sim_params->outfile);
	}
	fputc('\n',sim_params->outfile);
}

/* the N(i)-CA(i)-CA(i+1)-C(i+1) dihedral angles which constraints psi(i)+phi(i+1) together via the bias eta parameters (to determine their distributions) */
void cos_dihedral_naac(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	//TODO invalidate at chain breaks

	int i, ch;

	for (i = 1; i < chain->NAA-1; i++) {

		if (Distb(i, i) > 0. && Distb(i+1, i+1) > 0.){ /* alpha helix */
			ch = 'H';
		} else if (Distb(i, i) < 0. && Distb(i+1, i+1) < 0.){ /* beta strand */
			ch = 'E';
		} else {
			ch = 'C';
		}
		// naac dihedral angle's cosine
		vector x, y, z;
		subtract(x, chain->aa[i].ca, chain->aa[i].n);
		subtract(y, chain->aa[i+1].ca, chain->aa[i].ca);
		subtract(z, chain->aa[i+1].c, chain->aa[i+1].ca);
		fprintf(sim_params->outfile,"%c %g\n", ch, -cosdihedral(x, y, z));

	}
}

/* number of contacts between two amino acids */
void n_cont_aa( AA a, AA b, int *n_i_bb, int *n_i_sch, int *n_j_bb, int *n_j_sch, simulation_params *sim_params, int print_distances, int filter_out_14, int hbond_contact)
{

	/* exclude all 1-4 interactions */

	/* exclude interactions due to hydrogen bonding proximity */

	if (hbond_contact) return;

	*n_i_bb = 0;
	*n_i_sch = 0;
	*n_j_bb = 0;
	*n_j_sch = 0;

	//directed sequence distance, if on the same chain
	int seqdist;
	if (a.chainid != b.chainid)
		seqdist = 1000 * abs(b.chainid - a.chainid);
	else
		seqdist = b.num - a.num;

	char label1[4];
	char label2[4];
	if (a.etc & CA_) {
	    if (b.etc & CA_ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.ca,b.ca)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CA_--CA_ %g\n",sqrt(distance(a.ca,b.ca)));
		    *n_i_bb += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & N__ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.ca,b.n)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CA_--N__ %g\n",sqrt(distance(a.ca,b.n)));
		    *n_i_bb += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & C__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.ca,b.c)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CA_--C__ %g\n",sqrt(distance(a.ca,b.c)));
		    *n_i_bb += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & O__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.ca,b.o)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CA_--O_ %g\n",sqrt(distance(a.ca,b.o)));
		    *n_i_bb += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & CB_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.ca,b.cb)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CA_--CB_ %g\n",sqrt(distance(a.ca,b.cb)));
		    *n_i_bb += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.ca,b.g)) {
		if (print_distances) {
		    if (b.id=='C') {
			sprintf(label1,"S__");
		    } else if (b.id=='S' || b.id=='T') {
			sprintf(label1,"O__");
		    } else {
			sprintf(label1,"CB_");
		    }
		    fprintf(sim_params->outfile,"dist CA_--%3s_ %g\n",label1,sqrt(distance(a.ca,b.g)));
		}
		    *n_i_bb += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G2_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.ca,b.g2)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CA_--CB_ %g\n",sqrt(distance(a.ca,b.g2)));
		    *n_i_bb += 1;
		    *n_j_sch += 1;
		}
	    }
	}

	/*a->N*/
	if (a.etc & N__) {
	    if (b.etc & CA_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.n,b.ca)) {
		if (print_distances) fprintf(sim_params->outfile,"dist N__--CA_ %g\n",sqrt(distance(a.n,b.ca)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & N__ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.n,b.n)) {
		if (print_distances) fprintf(sim_params->outfile,"dist N__--N__ %g\n",sqrt(distance(a.n,b.n)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & C__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.n,b.c)) {
		if (print_distances) fprintf(sim_params->outfile,"dist N__--C__ %g\n",sqrt(distance(a.n,b.c)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & O__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.n,b.o)) {
		if (print_distances) fprintf(sim_params->outfile,"dist N__--O__ %g\n",sqrt(distance(a.n,b.o)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & CB_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
//		fprintf(stderr,"comparing %g > %g ?\n",sim_params->protein_model.touch2,distance(a.n,b.cb));
		if (sim_params->protein_model.touch2 > distance(a.n,b.cb)) {
		if (print_distances) fprintf(sim_params->outfile,"dist N__--CB_ %g\n",sqrt(distance(a.n,b.cb)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.n,b.g)) {
		if (print_distances) {
		    if (b.id=='C') {
			sprintf(label1,"S__");
		    } else if (b.id=='S' || b.id=='T') {
			sprintf(label1,"O__");
		    } else {
			sprintf(label1,"CB_");
		    }
		    fprintf(sim_params->outfile,"dist N__--%3s_ %g\n",label1,sqrt(distance(a.n,b.g)));
		}
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G2_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.n,b.g2)) {
		if (print_distances) fprintf(sim_params->outfile,"dist N__--CB_ %g\n",sqrt(distance(a.n,b.g2)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	}

	/*a->C*/
	if (a.etc & C__) {
	    if (b.etc & CA_ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.c,b.ca)) {
		if (print_distances) fprintf(sim_params->outfile,"dist C__--CA_ %g\n",sqrt(distance(a.c,b.ca)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & N__ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.c,b.n)) {
		if (print_distances) fprintf(sim_params->outfile,"dist C__--N__ %g\n",sqrt(distance(a.c,b.n)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & C__ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.c,b.c)) {
		if (print_distances) fprintf(sim_params->outfile,"dist C__--C__ %g\n",sqrt(distance(a.c,b.c)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & O__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.c,b.o)) {
		if (print_distances) fprintf(sim_params->outfile,"dist C__--O__ %g\n",sqrt(distance(a.c,b.o)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & CB_ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
//		fprintf(stderr,"comparing %g > %g ?\n",sim_params->protein_model.touch2,distance(a.c,b.cb));
		if (sim_params->protein_model.touch2 > distance(a.c,b.cb)) {
		if (print_distances) fprintf(sim_params->outfile,"dist C__--CB_ %g\n",sqrt(distance(a.c,b.cb)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.c,b.g)) {
		if (print_distances) {
		    if (b.id=='C') {
			sprintf(label1,"S__");
		    } else if (b.id=='S' || b.id=='T') {
			sprintf(label1,"O__");
		    } else {
			sprintf(label1,"CB_");
		    }
		    fprintf(sim_params->outfile,"dist C__--%3s_ %g\n",label1,sqrt(distance(a.c,b.g)));
		}
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G2_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.c,b.g2)) {
		if (print_distances) fprintf(sim_params->outfile,"dist C__--G2_ %g\n",sqrt(distance(a.c,b.g2)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	}

	/*a->O*/
	if (a.etc & O__) {
	    if (b.etc & CA_ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.o,b.ca)) {
		if (print_distances) fprintf(sim_params->outfile,"dist O__--CA_ %g\n",sqrt(distance(a.o,b.ca)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & N__ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.o,b.n)) {
		if (print_distances) fprintf(sim_params->outfile,"dist O__--N__ %g\n",sqrt(distance(a.o,b.n)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & C__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.o,b.c)) {
		if (print_distances) {
		    if ((b.num-a.num)==1) {
			fprintf(sim_params->outfile,"dist O__--C__ %g Y\n",sqrt(distance(a.o,b.c)));
		    } else {
			fprintf(sim_params->outfile,"dist O__--C__ %g\n",sqrt(distance(a.o,b.c)));
		    }
		}
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & O__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.o,b.o)) {
		if (print_distances) fprintf(sim_params->outfile,"dist O__--O__ %g\n",sqrt(distance(a.o,b.o)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & CB_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
//		fprintf(stderr,"comparing %g > %g ?\n",sim_params->protein_model.touch2,distance(a.o,b.cb));
		if (sim_params->protein_model.touch2 > distance(a.o,b.cb)) {
		if (print_distances) {
		    if ((b.num-a.num)==1) {
			fprintf(sim_params->outfile,"dist O__--CB_ %g X\n",sqrt(distance(a.o,b.cb)));
		    } else {
			fprintf(sim_params->outfile,"dist O__--CB_ %g\n",sqrt(distance(a.o,b.cb)));
		    }
		}
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G__) {
		if (sim_params->protein_model.touch2 > distance(a.o,b.g)) {
		if (print_distances) {
		    if (b.id=='C') {
			sprintf(label1,"S__");
		    } else if (b.id=='S' || b.id=='T') {
			sprintf(label1,"O__");
		    } else {
			sprintf(label1,"CB_");
		    }
		    fprintf(sim_params->outfile,"dist O__--%3s_ %g\n",label1,sqrt(distance(a.o,b.g)));
		}
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G2_) {
		if (sim_params->protein_model.touch2 > distance(a.o,b.g2)) {
		if (print_distances) fprintf(sim_params->outfile,"dist O__--CB_ %g\n",sqrt(distance(a.o,b.g2)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	}

	/*a->CB*/
	if (a.etc & CB_) {
	    if (b.etc & CA_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.cb,b.ca)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--CA_ %g\n",sqrt(distance(a.cb,b.ca)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & N__ && ( !filter_out_14 || (filter_out_14 && seqdist!=1 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.cb,b.n)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--N__ %g\n",sqrt(distance(a.cb,b.n)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & C__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.cb,b.c)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--C__ %g\n",sqrt(distance(a.cb,b.c)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & O__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.cb,b.o)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--O__ %g\n",sqrt(distance(a.cb,b.o)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & CB_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
//		fprintf(stderr,"comparing %g > %g ?\n",sim_params->protein_model.touch2,distance(a.cb,b.cb));
		if (sim_params->protein_model.touch2 > distance(a.cb,b.cb)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--CB_ %g\n",sqrt(distance(a.cb,b.cb)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.cb,b.g)) {
		if (print_distances) {
		    if (b.id=='C') {
			sprintf(label1,"S__");
		    } else if (b.id=='S' || b.id=='T') {
			sprintf(label1,"O__");
		    } else {
			sprintf(label1,"CB_");
		    }
		    fprintf(sim_params->outfile,"dist CB_--%3s_ %g\n",label1,sqrt(distance(a.cb,b.g)));
		}
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G2_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.cb,b.g2)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--CB_ %g\n",sqrt(distance(a.cb,b.g2)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	}

	/*a->G*/
	if (a.etc & G__) {
		if (print_distances) {
		    if (a.id=='C') {
			sprintf(label1,"S__");
		    } else if (a.id=='S' || a.id=='T') {
			sprintf(label1,"O__");
		    } else {
			sprintf(label1,"CB_");
		    }
		}
	    if (b.etc & CA_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g,b.ca)) {
		if (print_distances) fprintf(sim_params->outfile,"dist %3s--CA_ %g\n",label1,sqrt(distance(a.g,b.ca)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & N__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g,b.n)) {
		if (print_distances) fprintf(sim_params->outfile,"dist %3s--N__ %g\n",label1,sqrt(distance(a.g,b.n)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & C__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g,b.c)) {
		if (print_distances) fprintf(sim_params->outfile,"dist %3s--C__ %g\n",label1,sqrt(distance(a.g,b.c)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & O__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) { //we've already included it once
		if (sim_params->protein_model.touch2 > distance(a.g,b.o)) {
		if (print_distances) fprintf(sim_params->outfile,"dist %3s--O__ %g\n",label1,sqrt(distance(a.g,b.o)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & CB_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
//		fprintf(stderr,"comparing %g > %g ?\n",sim_params->protein_model.touch2,distance(a.g,b.cb));
		if (sim_params->protein_model.touch2 > distance(a.g,b.cb)) {
		if (print_distances) fprintf(sim_params->outfile,"dist %3s--CB_ %g\n",label1,sqrt(distance(a.g,b.cb)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g,b.g)) {
		if (print_distances) {
		    if (b.id=='C') {
			sprintf(label2,"S__");
		    } else if (b.id=='S' || b.id=='T') {
			sprintf(label2,"O__");
		    } else {
			sprintf(label2,"CB_");
		    }
		    fprintf(sim_params->outfile,"dist %3s--%3s_ %g\n",label1,label2,sqrt(distance(a.cb,b.g)));
		}
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G2_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g,b.g2)) {
		if (print_distances) fprintf(sim_params->outfile,"dist %3s--G2_ %g\n",label1,sqrt(distance(a.g,b.g2)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	}

	/*a->G2*/
	if (a.etc & G2_) {
	    if (b.etc & G2_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g2,b.ca)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--CA_ %g\n",sqrt(distance(a.g2,b.ca)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & N__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g2,b.n)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--N__ %g\n",sqrt(distance(a.g2,b.n)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & C__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g2,b.c)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--C__ %g\n",sqrt(distance(a.g2,b.c)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & O__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) { //we've already included it once
		if (sim_params->protein_model.touch2 > distance(a.g2,b.o)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--O__ %g\n",sqrt(distance(a.g2,b.o)));
		    *n_i_sch += 1;
		    *n_j_bb += 1;
		}
	    }
	    if (b.etc & CB_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
//		fprintf(stderr,"comparing %g > %g ?\n",sim_params->protein_model.touch2,distance(a.g2,b.cb));
		if (sim_params->protein_model.touch2 > distance(a.g2,b.cb)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--CB_ %g\n",sqrt(distance(a.g2,b.cb)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G__ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g2,b.g)) {
		if (print_distances) {
		    if (b.id=='C') {
			sprintf(label1,"S__");
		    } else if (b.id=='S' || b.id=='T') {
			sprintf(label1,"O__");
		    } else {
			sprintf(label1,"CB_");
		    }
		    fprintf(sim_params->outfile,"dist CB_--%3s_ %g\n",label1,sqrt(distance(a.g2,b.g)));
		}
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	    if (b.etc & G2_ && ( !filter_out_14 || (filter_out_14 && a.num!=b.num) )) {
		if (sim_params->protein_model.touch2 > distance(a.g2,b.g2)) {
		if (print_distances) fprintf(sim_params->outfile,"dist CB_--CB_ %g\n",sqrt(distance(a.g2,b.g2)));
		    *n_i_sch += 1;
		    *n_j_sch += 1;
		}
	    }
	}


//	fprintf(stderr,"calculated %d %d %d %d\n",*n_i_bb,*n_i_sch,*n_j_bb,*n_j_sch);

}

/* TODO: not just distance, but also different vdW radius for different atomic types and side chain atoms */
/* number of contacts per residue */
void number_of_contacts(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{

	fprintf(stderr,"Calculating number of contacts with touching distance %g\n",sqrt(sim_params->protein_model.touch2));

	int i, j;
	int *n_contacts_bb;
	int *n_contacts_sch;
	int n_i_bb, n_i_sch;
	int n_j_bb, n_j_sch;

	n_contacts_sch = (int *)calloc(chain->NAA,sizeof(int));
	n_contacts_bb = (int *)calloc(chain->NAA,sizeof(int));


	for (i = 1; i < chain->NAA-1; i++) {
	    for (j = i+1; j < chain->NAA-1; j++) {

		/* skip neighbouring amino acids on the same chain */
		if (j-i == 1 && chain->aa[i].chainid == chain->aa[j].chainid)
			continue;

		/* find all contacts between the two residues */
//		fprintf(stderr,"etc %x and %x\n",chain->aa[i].etc,chain->aa[j].etc);
		n_cont_aa(chain->aa[i],chain->aa[j],&n_i_bb,&n_i_sch,&n_j_bb,&n_j_sch, sim_params, /* print_distances = */ 0, /* filter_out_14 */ 0, 0);

		/* backbone contacts */
		n_contacts_bb[i] += n_i_bb;
		n_contacts_bb[j] += n_j_bb;
		/* side chain contacts */
		n_contacts_sch[i] += n_i_sch;
		n_contacts_sch[j] += n_j_sch;
	    }
	}

	/* print number of contacts with other details */
	fprintf(sim_params->outfile,"#n_cont_bb n_cont_sch resname(num) distb(i,i)\n");
	for (i = 1; i < chain->NAA-1; i++) {
	    fprintf(sim_params->outfile,"%d %d %d %g\n",n_contacts_bb[i],n_contacts_sch[i], chain->aa[i].id - 'A', Distb(i,i));
	}

}

/* print all the atomic distances for distribution generation */
void atomic_distances(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{

	fprintf(stderr,"Calculating atomic distances from i,i+2 (this excludes more than the 1-4 exclusion)!\n");

	int i, j;
	int n_i_bb, n_i_sch;
	int n_j_bb, n_j_sch;


	for (i = 1; i < chain->NAA-1; i++) {
//	    for (j = i; j < chain->NAA-1; j++) {
	    for (j = i+1; j < chain->NAA-1; j++) {

		/* skip neighbouring amino acids on the same chain */
		if (j-i == 1 && chain->aa[i].chainid == chain->aa[j].chainid)
			continue;

		/* calc if the residues are in contact due to a H-bond */
		/* to exclude the LJ interactions due to H-bond proximity */
		int hbond_proximity = 0;
		for (int i1 = i-1; i1<i+2; i1++) {
		    for (int j1 = j-1; j1<j+2; j1++) {
			if (((i1>=1) || (i1<chain->NAA)) && ((j1>=1) || (j1<chain->NAA))) {
			    if (( hstrength(chain->aa[i1].n,chain->aa[i1].h,chain->aa[j1].o,chain->aa[j1].c, &(sim_params->protein_model)) != 0 ) ||
				( hstrength(chain->aa[j1].n,chain->aa[j1].h,chain->aa[i1].o,chain->aa[i1].c, &(sim_params->protein_model)) != 0 )) {
				hbond_proximity = 1;
			    }
			}
		    }
		}

		/* find all contacts between the two residues */
		n_cont_aa(chain->aa[i],chain->aa[j],&n_i_bb,&n_i_sch,&n_j_bb,&n_j_sch, sim_params, /* print_distances = */ 1, /* filter_out_14 */ 1, hbond_proximity);
//		n_cont_aa(chain->aa[i],chain->aa[j],&n_i_bb,&n_i_sch,&n_j_bb,&n_j_sch, sim_params, /* print_distances = */ 1, /* filter_out_14 */ 1, 0);

	    }
	}

}

/* the CA-CA and CB-CB distances of residues in beta-beta contacts (to determine their distributions) */
//TODO: multi-chain proteins
void strand_bias_distances(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, j;

//fprintf(sim_params->outfile,"working hard on beta-beta contact.\n");
	/* loop over the bias map offdiagonal elements and find beta-beta contacts */
	for (i = 1; i < chain->NAA-1; i++) {
	    for (j = i+1; j < chain->NAA; j++) {

		if (Distb(i, i) < 0. && Distb(j, j) < 0. && Distb(i, j) > 0.) { /* beta-beta contact */
//fprintf(sim_params->outfile,"found beta-beta contact.\n");
		    if ((chain->aa[i].etc & CB_) && (chain->aa[j].etc & CB_)) {
			fprintf(sim_params->outfile,"E CBCB %g\n",sqrt(distance(chain->aa[i].cb,chain->aa[j].cb)));
		    }
		    if ((chain->aa[i].etc & CA_) && (chain->aa[j].etc & CA_)) {
			fprintf(sim_params->outfile,"E CACA %g\n",sqrt(distance(chain->aa[i].ca,chain->aa[j].ca)));
		    }
		}

	    }
	}

	/* loop over the bias map diagonal elements and find alpha contacts */
	for (i = 1; i < chain->NAA-3; i++) {

		if (Distb(i, i) > 0. && Distb(i+1, i+1) > 0. && Distb(i+2, i+2) > 0. && Distb(i+3, i+3) > 0. && Distb(i,i+3) > 0.) { /* alpha-alpha contact */
		    if ((chain->aa[i].etc & CB_) && (chain->aa[i+3].etc & CB_)) {
			fprintf(sim_params->outfile,"H CBCB %g\n",sqrt(distance(chain->aa[i].cb,chain->aa[i+3].cb)));
		    }
		    if ((chain->aa[i].etc & CA_) && (chain->aa[i+3].etc & CA_)) {
			fprintf(sim_params->outfile,"H CACA %g\n",sqrt(distance(chain->aa[i].ca,chain->aa[i+3].ca)));
		    }
		}

	}

}

/* the CA-CA and CB-CB distances of residues in beta-beta contacts (to determine their distributions) */
void hydrophobic_distances(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, j;
	char label1;
	char label2;

//	fprintf(sim_params->outfile,"%c %x\n",'A',hydrophobic_atoms_list('A',sim_params->protein_model.sidechain_properties));     //  CB_,             
//	fprintf(sim_params->outfile,"%c %x\n",'B',hydrophobic_atoms_list('B',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'C',hydrophobic_atoms_list('C',sim_params->protein_model.sidechain_properties));     //  CB_ & G__,       
//	fprintf(sim_params->outfile,"%c %x\n",'D',hydrophobic_atoms_list('D',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'E',hydrophobic_atoms_list('E',sim_params->protein_model.sidechain_properties));     //  G__,                                
//	fprintf(sim_params->outfile,"%c %x\n",'F',hydrophobic_atoms_list('F',sim_params->protein_model.sidechain_properties));     //  CB_ & G__,       
//	fprintf(sim_params->outfile,"%c %x\n",'G',hydrophobic_atoms_list('G',sim_params->protein_model.sidechain_properties));     //  0x0,             
//	fprintf(sim_params->outfile,"%c %x\n",'H',hydrophobic_atoms_list('H',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'I',hydrophobic_atoms_list('I',sim_params->protein_model.sidechain_properties));     //  CB_ & G__ & G2_, 
//	fprintf(sim_params->outfile,"%c %x\n",'K',hydrophobic_atoms_list('K',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'L',hydrophobic_atoms_list('L',sim_params->protein_model.sidechain_properties));     //  CB_ & G__,       
//	fprintf(sim_params->outfile,"%c %x\n",'M',hydrophobic_atoms_list('M',sim_params->protein_model.sidechain_properties));     //  CB_ & G__,       
//	fprintf(sim_params->outfile,"%c %x\n",'N',hydrophobic_atoms_list('N',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'P',hydrophobic_atoms_list('P',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'Q',hydrophobic_atoms_list('Q',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'R',hydrophobic_atoms_list('R',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'S',hydrophobic_atoms_list('S',sim_params->protein_model.sidechain_properties));     //  G__,             
//	fprintf(sim_params->outfile,"%c %x\n",'T',hydrophobic_atoms_list('T',sim_params->protein_model.sidechain_properties));     //  G2_,
//	fprintf(sim_params->outfile,"%c %x\n",'V',hydrophobic_atoms_list('V',sim_params->protein_model.sidechain_properties));     //  CB_ & G__ & G2_, 
//	fprintf(sim_params->outfile,"%c %x\n",'W',hydrophobic_atoms_list('W',sim_params->protein_model.sidechain_properties));     //  CB_ & G__,       
//	fprintf(sim_params->outfile,"%c %x\n",'Y',hydrophobic_atoms_list('Y',sim_params->protein_model.sidechain_properties));     //  CB_ & G__,       
//	fprintf(sim_params->outfile,"%c %x\n",'Z',hydrophobic_atoms_list('Z',sim_params->protein_model.sidechain_properties));     //  G__,             
//
//	return;

	/* loop over the residues and find the hydrophobic (HH/HA/HP) contacts */
	for (i = 1; i < chain->NAA-1; i++) {
	    for (j = i+1; j < chain->NAA; j++) {

		/* set label */
		if ((chain->aa[i].etc & HYDROPHOBIC) && (chain->aa[i].etc & HYDROPHOBIC)) { /* H-H */
		   label1 = 'H';
		   label2 = 'H';
		} else if ( ((chain->aa[i].etc & HYDROPHOBIC) && (chain->aa[j].etc & AMPHIPATHIC)) ||
			    ((chain->aa[i].etc & AMPHIPATHIC) && (chain->aa[j].etc & HYDROPHOBIC)) ) { /* H-A */
		   label1 = 'H';
		   label2 = 'A';
		} else if ((chain->aa[i].etc & HYDROPHOBIC) || (chain->aa[j].etc & HYDROPHOBIC)) { /* H-P */
		   label1 = 'H';
		   label2 = 'P';
		} else {
		   continue;
		}

		/* reorder in case to have the first one as hydrophobic */
		int ii, jj;
		if (chain->aa[i].etc & HYDROPHOBIC) {
			ii = i;
			jj = j;
		} else { /* j was hydrophobic */
			ii = j;
			jj = i;
		}

		/* calculate atomic distances */
		if ((hydrophobic_atoms_list(chain->aa[ii].id,sim_params->protein_model.sidechain_properties) & CB_) && (chain->aa[ii].etc & CB_)) {
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & CB_) && (chain->aa[jj].etc & CB_)) {
			fprintf(sim_params->outfile,"%c%c %c %c CBCB %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].cb,chain->aa[jj].cb)));
		    }
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & G__) && (chain->aa[jj].etc & G__)) {
			fprintf(sim_params->outfile,"%c%c %c %c CBG1 %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].cb,chain->aa[jj].g)));
		    }
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & G2_) && (chain->aa[jj].etc & G2_)) {
			fprintf(sim_params->outfile,"%c%c %c %c CBG2 %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].cb,chain->aa[jj].g2)));
		    }
		}
		if ((hydrophobic_atoms_list(chain->aa[ii].id,sim_params->protein_model.sidechain_properties) & G__) && (chain->aa[ii].etc & G__)) {
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & CB_) && (chain->aa[jj].etc & CB_)) {
			fprintf(sim_params->outfile,"%c%c %c %c G1CB %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].g,chain->aa[jj].cb)));
		    }
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & G__) && (chain->aa[jj].etc & G__)) {
			fprintf(sim_params->outfile,"%c%c %c %c G1G1 %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].g,chain->aa[jj].g)));
		    }
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & G2_) && (chain->aa[jj].etc & G2_)) {
			fprintf(sim_params->outfile,"%c%c %c %c G1G2 %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].g,chain->aa[jj].g2)));
		    }
		}
		if ((hydrophobic_atoms_list(chain->aa[ii].id,sim_params->protein_model.sidechain_properties) & G2_) && (chain->aa[ii].etc & G2_)) {
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & CB_) && (chain->aa[jj].etc & CB_)) {
			fprintf(sim_params->outfile,"%c%c %c %c G2CB %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].g2,chain->aa[jj].cb)));
		    }
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & G__) && (chain->aa[jj].etc & G__)) {
			fprintf(sim_params->outfile,"%c%c %c %c G2G1 %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].g2,chain->aa[jj].g)));
		    }
		    if ((hydrophobic_atoms_list(chain->aa[jj].id,sim_params->protein_model.sidechain_properties) & G2_) && (chain->aa[jj].etc & G2_)) {
			fprintf(sim_params->outfile,"%c%c %c %c G2G2 %g\n",label1,label2,chain->aa[ii].id,chain->aa[jj].id,sqrt(distance(chain->aa[ii].g2,chain->aa[jj].g2)));
		    }
		}

	    }
	}

}

/*calc the half angle between CB-CA-HA*/
double half_cb_ca_ha_angle(AA *a) {

	vector can, cac, cn, cacb;

	subtract(can,a->n,a->ca);
	subtract(cac,a->c,a->ca);
	normalize(can);
	normalize(cac);

	add(cn, can, cac);

	subtract(cacb,a->cb,a->ca);

	return angle(cacb,cn);
}

/* Ramachandran conformation of all residues */
//TODO: multi-chain proteins
void phipsi(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i;
	double phi, psi, tau, chi1, chi2, cb_ca_angle;
	vector cac, can;
	char *omega;
	vector x,y,z;
	double dihedral_etaA_NAAC;
	double dihedral_etaB_NAAC;
	double dihedral_etaA_BAAB;
	double dihedral_etaB_BAAB;

	for (i = 1; i < chain->NAA; i++) {
		all_torsions(chain,biasmap,&phi, &psi, &chi1, &chi2, i);

		subtract(cac, chain->aa[i].c, chain->aa[i].ca);
		subtract(can, chain->aa[i].n, chain->aa[i].ca);
		tau = M_180_PI * angle(cac, can);

		omega = (chain->aa[i].etc & CIS) ? "  0" : "180";

		if (chain->aa[i].etc & CB_) {
			cb_ca_angle = M_180_PI * half_cb_ca_ha_angle(&(chain->aa[i]));
		} else {
			cb_ca_angle = 0./0.;
		}

		/* for beta or alpha residues print N(i)-CA(i)-CA(i+1)-C(i+1) dihedral */
		dihedral_etaA_NAAC = 0;
		dihedral_etaB_NAAC = 0;
		dihedral_etaA_BAAB = 0;
		dihedral_etaB_BAAB = 0;
		if (i<chain->NAA-1) {
		    /* N-CA-CA-C */
		    subtract(x, chain->aa[i].ca, chain->aa[i].n);
		    subtract(y, chain->aa[i+1].ca, chain->aa[i].ca);
		    subtract(z, chain->aa[i+1].c, chain->aa[i+1].ca);
		    if (Distb(i,i+1) < 0) { /* beta twist */
			dihedral_etaB_NAAC = -acos(cosdihedral(x, y, z)); 
		    } else if (Distb(i,i+1) > 0) { /* alpha twist */
			//dihedral_etaA_NAAC = -phasindihedral(x,y,z, 0.13917, 0.99); 
			dihedral_etaA_NAAC = -acos(cosdihedral(x, y, z)); 
		    }
		    /* CB-CA-CA-CB */
		    if (chain->aa[i].etc & CB_ && chain->aa[i+1].etc & CB_) {
			subtract(x, chain->aa[i].ca, chain->aa[i].cb);
			subtract(y, chain->aa[i+1].ca, chain->aa[i].ca);
			subtract(z, chain->aa[i+1].cb, chain->aa[i+1].ca);
			if (Distb(i,i+1) < 0) { /* beta twist */
			    dihedral_etaB_BAAB = -acos(cosdihedral(x, y, z)); 
			} else if (Distb(i,i+1) > 0) { /* alpha twist */
			    //dihedral_etaA_BAAB = -phasindihedral(x,y,z, 0.13917, 0.99); 
			    dihedral_etaA_BAAB = -acos(cosdihedral(x, y, z)); 
			}
		    }
		}

		fprintf(sim_params->outfile,"%c%7.1f%7.1f%7.1f  %s %7.1f%7.1f%7.1f%7.3f%7.3f%7.3f%7.3f\n", chain->aa[i].id, phi, psi, tau,
		       omega, chi1, chi2, cb_ca_angle, dihedral_etaA_NAAC, dihedral_etaB_NAAC, dihedral_etaA_BAAB, dihedral_etaB_BAAB);
	}
}

/* Print the geometrical data for every single hydrogen bond */
void hbond_geometry(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm) {

	int i, j;
	vector oh;
	vector hn;
	vector co;

	for (i = 1; i < chain->NAA; i++) {
	    for (j = 1; j < chain->NAA; j++) {

		AA *a = &(chain->aa[i]);
		AA *b = &(chain->aa[j]);
		if (!(a->etc & H__ && b->etc & O__)) continue;

		subtract(oh, a->h, b->o);
		if (square(oh) > sim_params->protein_model.hboh2) continue; // d(O,H) is too long
		subtract(hn, a->n, a->h);
		if (cosine(oh, hn) < sim_params->protein_model.hbohn) continue; // <(O,H,N) is too small
		subtract(co, b->o, b->c);
		if (cosine(co, oh) < sim_params->protein_model.hbcoh) continue; // <(C,O,H) is too small
		/* found a H-bond */
		fprintf(sim_params->outfile,"%d %d %g %g %g\n",i,j,sqrt(square(oh)),cosine(oh,hn),cosine(co,oh));
	    }
	}
}

/* Number of hydrogen bonds */
void hbtot(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, j, tot = 0;

	for (i = 1; i < chain->NAA; i++)
		for (j = 1; j < chain->NAA; j++)
			tot += hdonor((chain->aa) + i, (chain->aa) + j, &(sim_params->protein_model));

	fprintf(sim_params->outfile,"Hbonds = %d\n", tot);
}

/* Hydrogen bonding patterns (within the same chain) */
void hpattern(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, j;
	unsigned int donor, accep;
	AA *a, *b;

	int *forward_pattern;
	int *backward_pattern;

	forward_pattern = (int *)malloc(chain->NAA * sizeof(int));
	backward_pattern = (int *)malloc(chain->NAA * sizeof(int));
	for (i = 0; i< chain->NAA; i++) {
		forward_pattern[i] = 0;
		backward_pattern[i] = 0;
	}


	for (i = 5; i < chain->NAA - 4; i++) {

		donor = accep = 0;
		a = chain->aa + i;

		for (j = 1; j < chain->NAA; j++) {

			if (i == j)
				continue;
			b = (chain->aa) + j;

			/* H-pattern within the same chain only */
			if (a->chainid != b->chainid)
				continue;


			if (hdonor(a, b, &(sim_params->protein_model))) {
				//fprintf(stderr,"Found hbond don: %d to acc: %d\n",i,j);
				donor |= 1 << (15 - j + i);
				if (i<j) forward_pattern[j-i-1]++;
				else backward_pattern[i-j-1]++;
			}
			if (hdonor(b, a, &(sim_params->protein_model))) {
				//fprintf(stderr,"Found hbond don: %d to acc: %d\n",j,i);
				accep |= 1 << (15 - i + j);
				if (i<j) backward_pattern[j-i-1]++;
				else forward_pattern[i-j-1]++;
			}
		}
		//fprintf(sim_params->outfile,"%3d   %.8X   %.8X\n", i, donor, accep);
		//fprintf(sim_params->outfile,"%3d   %8X   %8X\n", i, donor, accep);
	}
	fprintf(sim_params->outfile,"Forward_pattern ");
	for (i = 0; i< chain->NAA; i++) {
		fprintf(sim_params->outfile,"%d ", forward_pattern[i]);
	}
	fprintf(sim_params->outfile,"\n");
	fprintf(sim_params->outfile,"Backward_pattern ");
	for (i = 0; i< chain->NAA; i++) {
			fprintf(sim_params->outfile,"%d ", backward_pattern[i]);
		}
	fprintf(sim_params->outfile,"\n");

	if (forward_pattern) free(forward_pattern);
	if (backward_pattern) free(backward_pattern);
}

void stepsize(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i;
	double phi, psi, dd, size = 0.0;
	static double ph[128], ps[128];

	//static int swit = -1;

	for (i = 1; i < chain->NAA && i < 128; i++) {
		torsion(chain,biasmap,&phi, &psi, i);
        //if(i>1){
		dd = phi - ph[i];
		ph[i] = phi;
		if (dd > 180.0)
			dd -= 360.0;
		if (dd < -180.0)
			dd += 360.0;
		size += dd * dd;
	   //}
		//if(i < NAA - 1){
		
		dd = psi - ps[i];
		ps[i] = psi;
		if (dd > 180.0)
			dd -= 360.0;
		if (dd < -180.0)
			dd += 360.0;
		size += dd * dd;
	     //}
	}

   //if(swit == 1){
    //fprintf(sim_params->outfile,"Amplitude = %g   Acceptance = %g   Drift = %g ",  amplitude, acceptance, sqrt(size));
    //}
    //swit *= -1;
	
  fprintf(sim_params->outfile,"Amplitude = %g   Acceptance = %g   Drift = %g\n",  sim_params->amplitude, sim_params->acceptance, sqrt(size));
	
}

void energy_gradient_wrt_parameters(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	energy_probe_1(chain,biasmap,sim_params);
//#ifndef INTERNAL_CDLEARN
	/* print all the gradients */
	for (int i=0; i<34; i++) {
	    fprintf(sim_params->outfile,"%f ", sim_params->energy_gradient[i]);
	}
	fputc('\n',sim_params->outfile);
//#endif

}

/* number of contacts */
int count_contacts(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, int native)
{
	int i, j, contax;
	if (native && biasmap->distb == NULL)
		return -1;

	contax = 0;
	for (i = 1; i < chain->NAA; i++) {
		for (j = 1; j < i; j++) {

			/* skip neighbours on the same chain */
			if (j==i-1 && chain->aa[i].chainid == chain->aa[j].chainid)
				continue;

			if (contact((chain->aa) + i, (chain->aa) + j, &(sim_params->protein_model))) {
				if (native) {
					if (Distb(i, j) > 0.0)
						contax++;
				} else {
					contax++;
				}
			}
		}
	}

	return contax;
}

/* number of native contacts */
void n_native_contacts(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int n_cont = count_contacts(chain,biasmap,sim_params, /* native = */ 1);
	if (n_cont < 0) {
	    return;
	} else {
	    fprintf(sim_params->outfile,"Number of native contacts = %i\n", n_cont);
	}

}

/* number of all contacts */
void n_contacts(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int n_cont = count_contacts(chain,biasmap,sim_params, /* native = */ 0);
	fprintf(sim_params->outfile,"Number of contacts = %i\n", n_cont);

}

void pdbout(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	//fprintf(stderr,"Outputting PDB, %d amino acids.\n", chain->NAA-1);
	double E_tot = totenergy(chain);
	pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), sim_params->outfile, &E_tot);
/*  vector mol_com;
  mol_com[0] = mol_com[1] = mol_com[2] = 0.0;
  for (int i = 1; i < chain->NAA; i++){
	vector com;
	add(com, ((chain->aa) + i)->ca, ((chain->aa) + i)->n);
	add(com, com, ((chain->aa) + i)->c);
	scale(com,1.0/3.0, com);
	add(mol_com, mol_com, com);
  }
  scale(mol_com, 1.0/(double)(chain->NAA-1), mol_com);
  fprintf(stderr,"Centre-of-mass: %g %g %g\n",mol_com[0],mol_com[1],mol_com[2]); */
}


/* Print a snapshot in CHK format */
void checkpoint_out(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void*mpi_comm) {
	int current_stored = 1;

    if(sim_params->num_NS_per_checkpoint == 0) return;

    fprintf(stderr," iter %d num_NS %d ",sim_params->iter,sim_params->num_NS_per_checkpoint);

//    if( sim_params->iter % sim_params->num_NS_per_checkpoint == 0 ){
	update_sim_params_from_chain(chain,sim_params); // updating NAA and seq
	fprintf(stderr,"printing checkpoint file as test.\n");
	output_checkpoint_file(chain, current_stored, sim_params,mpi_comm);
//	}

}


/* contact map in text format, which is pretty wasteful */
/* WARNING!  The output of this contact map cannot be used in simulations with ./peptide.
	Fscanf skips space characters, and complains about lopsided bias.
	Even so, the amino acid names in the diagonal would be used to read in the secondary
	structure, according to the last 2 bits of the UTF-8 symbols, giving meaningless
	bias potential energy contributions. */
/* cm_ideal should be used instead. */
//TODO: multi-chain proteins
void cm_txt(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, j, ch;

	for (i = 1; i < chain->NAA; i++) {
		for (j = 1; j < chain->NAA; j++) {
			switch (abs(i - j)) {
			case 0:
				ch = chain->aa[j].id;
				break;
			case 1:
				ch = aligned((chain->aa) + i, (chain->aa) + j) ? '+' : '-';
				break;
			default:
				ch = contact((chain->aa) + i, (chain->aa) + j, &(sim_params->protein_model)) ? '#' : ' ';
			}
			fputc(ch,sim_params->outfile);
		}
		fputc('\n',sim_params->outfile);
	}
}

/* contact map as an image in PBM format */
void cm_pbm(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, j, byte, bit;

	fprintf(sim_params->outfile,"P4 %u %u\n", chain->NAA - 1, chain->NAA - 1);

	for (i = 1; i < chain->NAA; i++)
		for (j = 1; j < chain->NAA;) {
			for (byte = 0, bit = 0x80; bit; bit >>= 1, j++)
				if (j < chain->NAA && contact((chain->aa) + i, (chain->aa) + j, &(sim_params->protein_model)))
					byte |= bit;
			fputc(byte,sim_params->outfile);
		}

	fputc('\n',sim_params->outfile);
}

static int dumps(char *s,FILE *outfile)
{
	int retv;

	retv = fputs(s,outfile);
	fputc('\n',outfile);
	for (; *s; s++)
		*s = ' ';

	return retv;
}

//TODO: multi-chain protein
static void hbss_loose(Chain *chain,Biasmap *biasmap, int alpha_i_i4_contacts, simulation_params *sim_params)
{
	int i, j, ch, helix = 0;
	char *ss = NULL;
	char *map = NULL;

	ss = (char *) realloc(ss, (chain->NAA + 1) * sizeof(char));
	if (ss == NULL)
		return;

	ss[chain->NAA] = '\0';		/* terminate the string */
	/* allocate the whole contact map */
	map = (char *) realloc(map, (chain->NAA + 1) * (chain->NAA + 1) * sizeof(char));
	if (map == NULL)
		return;


	/* detect (3- and 4-) helices and coil segments */
	for (i = 1; i < chain->NAA; i++) {
		for (j = i + 3; j <= i + 4; j++)
			if (j < chain->NAA)
				helix += Hbm(j - 1, i);

		/* important to have this in-between */
		ss[i] = (helix > 1) ? 'H' : 'C';

		for (j = i - 4; j <= i - 3; j++)
			if (j > 0)
				helix -= Hbm(i - 1, j);
	}

	/* map beta-bridges and assign beta-stands */
	for (i = 1; i < chain->NAA; i++) {
		for (j = 1; j < chain->NAA; j++) {
			if (i == j) {
				ch = chain->aa[i].id;
				Map(i,j) = ch;
			} else if (Hbm(i - 1, j - 1) && Hbm(i, j)) {
				ch = '<';	/* helix */
				Map(i,j) = 'X';
			} else if (Hbm(j - 1, i - 1) && Hbm(j, i)) {
				ch = '>';	/* helix */
				Map(i,j) = 'X';
			} else if ((Hbm(i - 1, j - 1) && Hbm(j, i)) ||
				 (Hbm(j - 1, i - 1) && Hbm(i, j))) {
				ch = '\\';
				Map(i,j) = 'X';
				ss[i] = ss[j] = 'E';
			} else if ((Hbm(i, j - 1) && Hbm(j, i - 1)) ||
				   (Hbm(j - 1, i) && Hbm(i - 1, j))) {
				ch = '/';
				Map(i,j) = 'X';
				ss[i] = ss[j] = 'E';
			}
			else {
				ch = ' ';
				Map(i,j) = 'O';
			}
			//fputc(ch,sim_params->outfile);
		}
		//fputc('\n',sim_params->outfile);
	}
	//fputs(ss + 1,sim_params->outfile);

	/* Fill in diagonal and i,i+/-1 of the contact map
	   and print the contact map */
	/* Optionally, add alpha i,i+4 contacts */
	for (i = 1; i < chain->NAA; i++) {
		for (j = 1; j < chain->NAA; j++) {
			/* mark secondary structure in diagonal */
			if (i == j) {
				switch (ss[i]) {
				case 'H': /* helix */
					Map(i,j) = 'U';
					break;
				case 'E': /* strand */
					Map(i,j) = 'Z';
					break;
				case 'C': /* coil */
					Map(i,j) = 'O';
				}
			/* mark secondary structure in i,i+1 */
			} else if (abs(i-j)==1 && ss[i] == ss[j]) {
				switch (ss[i]) {
				case 'H': /* helix */
					Map(i,j) = 'U';
					break;
				case 'E': /* strand */
					Map(i,j) = 'Z';
					break;
				case 'C': /* coil */
					Map(i,j) = 'O';
				}
			/* optionally mark alpha-helix in i,i+4 */
			} else if ( abs(i-j)==4 && alpha_i_i4_contacts && ss[i] == 'H' ) {
				if (i<j) { /* upper triangle */
					if ( Map(i,j-1) == 'X' && Map(i+1,j) == 'X' ) Map(i,j) = 'X';
				} else { /* lower triangle */
					if ( Map(i-1,j) == 'X' && Map(i,j+1) == 'X' ) Map(i,j) = 'X';
				}
			}
			fputc(Map(i,j),sim_params->outfile);
		}
		fputc('\n',sim_params->outfile);
	}
	free(ss);
	free(map);
}

//TODO: multi-chain protein
static void hbss_compact(Chain *chain,Biasmap *biasmap, simulation_params *sim_params)
{
	int i, j, k, flag;
	char *ss = NULL;
	char *sss;

	ss = (char *) realloc(ss, (2 * chain->NAA + 1) * sizeof(char));
	if (ss == NULL)
		return;

	fprintf(sim_params->outfile,"HBSS Report, N = %d\n", chain->NAA - 1);

	sss = ss + chain->NAA;
	sss[chain->NAA] = ss[chain->NAA] = '\0';	/* terminate the string */

	for (i = 1; i < chain->NAA; i++) {
		sss[i] = chain->aa[i].id;
		ss[i] = 'C';
	}
	dumps(sss + 1,sim_params->outfile);

	/* map helical segments */
	for (k = 2; k < chain->NAA - 1; k++) {
		flag = 0;
		for (i = 1, j = i + k; j < chain->NAA; i++, j++)
			if (Hbm(j - 1, i - 1) && Hbm(j, i)) {
				sss[i] = sss[i] == ' ' ? '>' : 'X';
				sss[j] = sss[j] == ' ' ? '<' : 'X';
				ss[i] = ss[j] = 'H';
				flag = 1;
			} else if (Hbm(i - 1, j - 1) && Hbm(i, j)) {
				sss[i] = sss[i] == ' ' ? '<' : 'O';
				sss[j] = sss[j] == ' ' ? '>' : 'O';
				ss[i] = ss[j] = 'H';
				flag = 1;
			} else if (flag) {
				dumps(sss + 1,sim_params->outfile);
				flag = 0;
			}

		if (flag)
			dumps(sss + 1,sim_params->outfile);
	}

	for (k = 5; k < 2 * chain->NAA - 4; k++) {
		flag = 0;
		i = (k - 3) >> 1;
		j = (k + 4) >> 1;
		for (; 0 < i && j < chain->NAA; i--, j++)
			if (Hbm(j - 1, i) && Hbm(j, i - 1)) {
				sss[i] = sss[i] == ' ' ? ')' : 'X';
				sss[j] = sss[j] == ' ' ? '(' : 'X';
				ss[i] = ss[j] = 'H';
				flag = 1;
			} else if (Hbm(i, j - 1) && Hbm(i - 1, j)) {
				sss[i] = sss[i] == ' ' ? '(' : 'O';
				sss[j] = sss[j] == ' ' ? ')' : 'O';
				ss[i] = ss[j] = 'H';
				flag = 1;
			} else if (flag) {
				dumps(sss + 1,sim_params->outfile);
				flag = 0;
			}

		if (flag)
			dumps(sss + 1,sim_params->outfile);
	}

	/* map parallel beta-sheets and helices */
	for (k = 2; k < chain->NAA - 1; k++) {
		flag = 0;
		for (i = 1, j = i + k; j < chain->NAA; i++, j++)
			if (Hbm(i - 1, j - 1) && Hbm(j, i)) {
				sss[i] = '}';
				sss[j] = ']';
				ss[i] = ss[j] = 'E';
				flag = 1;
			} else if (Hbm(j - 1, i - 1) && Hbm(i, j)) {
				sss[i] = ']';
				sss[j] = '}';
				ss[i] = ss[j] = 'E';
				flag = 1;
			} else if (flag) {
				dumps(sss + 1,sim_params->outfile);
				flag = 0;
			}

		if (flag)
			dumps(sss + 1,sim_params->outfile);
	}

	/* map anti-parallel beta-sheets */
	for (k = 5; k < 2 * chain->NAA - 4; k++) {
		flag = 0;
		i = (k - 3) >> 1;
		j = (k + 4) >> 1;
		for (; 0 < i && j < chain->NAA; i--, j++)
			if (Hbm(i, j - 1) && Hbm(j, i - 1)) {
				sss[i] = '[';
				sss[j] = ']';
				ss[i] = ss[j] = 'E';
				flag = 1;
			} else if (Hbm(j - 1, i) && Hbm(i - 1, j)) {
				sss[i] = '{';
				sss[j] = '}';
				ss[i] = ss[j] = 'E';
				flag = 1;
			} else if (flag) {
				dumps(sss + 1,sim_params->outfile);
				flag = 0;
			}

		if (flag)
			dumps(sss + 1,sim_params->outfile);
	}

	dumps(ss + 1,sim_params->outfile);
	free(ss);
}

/* Low level routine for DSSP-like H-bond-based secondary structure analysis */
/* secondary_structure=1:	output secondary structure
   contact_map=1:		output ideal contact map
   alpha_i_i4_contacts:		include i,i+4 alpha contacts (only for the contact map generation) */
//TODO: multi-chain protein
void dssp(Chain *chain,Biasmap *biasmap, int secondary_structure, int contact_map, int alpha_i_i4_contacts, simulation_params *sim_params) {
	int i, j;

	hbm = NULL;
	hbm = (int *) realloc(hbm, chain->NAA * chain->NAA * sizeof(int));
	if (hbm == NULL)
		return;

	/* blank out the map */
	for (i = 0; i < chain->NAA * chain->NAA; i++)
		hbm[i] = 0;

	/* map long-range hydrogen bonds between peptide bonds,
	   which can be joined by one and only one hydrogen bond,
	   unlike amino acids themselves */
	for (i = 1; i < chain->NAA; i++)
		for (j = 1; j < chain->NAA; j++)
			if (hdonor((chain->aa) + i, (chain->aa) + j, &(sim_params->protein_model)))
				Hbm(i - 1, j) = 1;

	/* for (i = 0; i < NAA; fputc('\n',sim_params->outfile), i++)
	   for (j = 0; j < NAA; j++)
	   fputc(Hbm(i, j) ? 'X' : ' ',sim_params->outfile); */

	/* tolerate a single missing bond in a pattern */
	for (i = 1; i < chain->NAA - 1; i++)
		for (j = 1; j < chain->NAA - 1; j++) {
			if (Hbm(i - 1, j - 1) && Hbm(i + 1, j + 1)) {
				if (abs(i - j) < 5)
					Hbm(i, j) = 1;	/* helix */
				else if (!Hbm(i, j))
					Hbm(j, i) = 1;	/* parallel */
			}
			if (Hbm(i - 1, j + 1) && Hbm(i + 1, j - 1)
			    && !Hbm(i, j))
				Hbm(j, i) = 1;	/* antiparallel */
		}

	if (secondary_structure) hbss_compact(chain,biasmap,sim_params);
	if (contact_map) hbss_loose(chain,biasmap,alpha_i_i4_contacts, sim_params);
	free(hbm);

}

/* DSSP-like H-bond-based secondary structure analysis */
void hbss(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm) {

	dssp(chain,biasmap,/* secondary_structure */ 1, /* contact_map */ 0, /* alpha_i_i4_contacts */ 0, sim_params);

}

/* DSSP-like H-bond-based ideal contact map generation without i,i+4 alpha contacts */
void cm_ideal(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm) {

	dssp(chain,biasmap,/* secondary_structure */ 0, /* contact_map */ 1, /* alpha_i_i4_contacts */ 0, sim_params);
	
}

/* DSSP-like H-bond-based ideal contact map generation with i,i+4 alpha contacts */
void cm_ideal_4(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm) {

	dssp(chain,biasmap,/* secondary_structure */ 0, /* contact_map */ 1, /* alpha_i_i4_contacts */ 1, sim_params);

}

void ergtot(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{	/*changed, can be removed */ 
	double tote = totenergy(chain);
	fprintf(sim_params->outfile, "Energy = totalE %.6f ( diagnolE %.6f extE %.6f firstlastE %.6f) Rotmers (%i)\n", tote, locenergy(chain), extenergy(chain), firstlastenergy(chain), bestRot(chain));
	fprintf(stderr, "totalE %.6f extE %.6f \n", tote, extenergy(chain));
}

void vdw_max_gamma(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm) {

//	char maxgamma_seq[27] = "ABCDEFGHIGKLMNGPQRSTGVWGYZ";
//	Chain *chain_maxgamma_test = NULL;
//	Chaint *chaint = NULL;
//
//	build_peptide_from_sequence(chain_maxgamma_test,chaint,maxgamma_seq, &sim_params);

	/* backbone constants */
//	vdw_backbone_constants(chain,&(sim_params->protein_model), sim_params->outfile, /* verbose = */ 1);
	/* gamma-gamma and gamma-nongamma constants */
//	vdw_maxgamma_calc(chain, &(sim_params->protein_model), sim_params->outfile, /* verbose = */ 0);
//	print_vdw_cutoff_distances(&(sim_params->protein_model),sim_params->outfile);

	vdw_cutoff_distances_calculate(sim_params, sim_params->outfile, 1);
}

void energy_contributions(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm){
  energy_contributions_in_energy_c(chain,biasmap,totenergy(chain), &(sim_params->protein_model), sim_params->outfile);
}

void exclude_energy_contributions(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm){
  exclude_energy_contributions_in_energy_c(chain,biasmap,totenergy(chain), &(sim_params->protein_model), sim_params->outfile);
}

void wrinkle(struct phasor *a, vector zi, vector p, vector q, vector zk)
{
	vector r, ri, rk;
	struct phasor b;

	subtract(r, q, p);

	crossprod(ri, zi, r);
	crossprod(rk, r, zk);

	b.y = dotprod(ri, zk) * sqrt(square(r));
	b.x = dotprod(ri, rk);
	b.k = 0;

	*a = phasiply(*a, b);
}

/* compact writhe calculator */
//TODO: multi-chain protein
void writhe(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm)
{
	int i, k;
	vector zi, zk;
	struct phasor a = { 0., 1., 0 };

	for (i = 1; i < chain->NAA - 3; i++) {
		subtract(zi, chain->aa[i + 1].ca, chain->aa[i].ca);
		for (k = i + 2; k < chain->NAA - 1; k++) {
			subtract(zk, chain->aa[k + 1].ca, chain->aa[k].ca);
			
			wrinkle(&a, zi, chain->aa[i].ca, chain->aa[k].ca, zk);
			wrinkle(&a, zi, chain->aa[k + 1].ca, chain->aa[i].ca, zk);
			wrinkle(&a, zi, chain->aa[i + 1].ca, chain->aa[k + 1].ca, zk);
			wrinkle(&a, zi, chain->aa[k].ca, chain->aa[i + 1].ca, zk);			
			rephase(&a);
	        
		}
	} 
	fprintf(sim_params->outfile,"Writhe = %g\n", phase(a) / (2 * M_PI)); 
}

/* explicit eigenvalues of symmetric 3x3 matrix given by two vectors,
implemented according to Oliver Smith 1961 */
static void eigensmith(vector eig, vector d, vector c)
{
	double m, p, q, phi;

	m = (d[0] + d[1] + d[2]) / 3;
	d[0] -= m;
	d[1] -= m;
	d[2] -= m;

	/* half-determinant */
	q = 0.5 * (d[0] * d[1] * d[2] - d[0] * c[0] * c[0] -
		   d[1] * c[1] * c[1] - d[2] * c[2] * c[2]) +
	    c[0] * c[1] * c[2];

	/* one sixth frobenius */
	p = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
	     + 2 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2])) / 6;

	/* Cardano's solution */
	phi = atan2(sqrt(p * p * p - q * q), q) / 3;

	p = sqrt(p);
	q = M_SQRT3 * p * sin(phi);
	p = p * cos(phi);

	eig[0] = m + 2 * p;
	eig[1] = m - p + q;
	eig[2] = m - p - q;
}

/* centroid and radius of gyration for alpha-carbon trace */
void rgyr(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm)
{
	int i;
	double w;
	vector ri, rc = { 0.0, 0.0, 0.0 };
	vector d = { 0.0, 0.0, 0.0 }, c = {
	0.0, 0.0, 0.0};

	w = 1.0 / (chain->NAA - 1);

	for (i = 1; i < chain->NAA; i++)
		add(rc, rc, chain->aa[i].ca);

	scale(rc, w, rc);

	for (i = 1; i < chain->NAA; i++) {
		subtract(ri, chain->aa[i].ca, rc);
		d[0] += ri[0] * ri[0];
		d[1] += ri[1] * ri[1];
		d[2] += ri[2] * ri[2];
		c[0] += ri[1] * ri[2];
		c[1] += ri[2] * ri[0];
		c[2] += ri[0] * ri[1];
	}

	scale(d, w, d);
	scale(c, w, c);

	w = sqrt(d[0] + d[1] + d[2]);

	eigensmith(ri, d, c);	/* this messes with d */

	fprintf(sim_params->outfile,"Rgyr = %8.3f   PrinComps = %8.3f%8.3f%8.3f      N = %4d\n",
	       w, sqrt(ri[0]), sqrt(ri[1]), sqrt(ri[2]), chain->NAA - 1);
	fprintf(sim_params->outfile,"Centroid      CA              %8.3f%8.3f%8.3f\n",
	       rc[0], rc[1], rc[2]);
}

/* centroid and radius of gyration for alpha-carbon trace */
void hydrophobic_rgyr(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm)
{
	int i;
	double w;
	vector ri, rc = { 0.0, 0.0, 0.0 };
	vector d = { 0.0, 0.0, 0.0 }, c = {
	0.0, 0.0, 0.0};

	int *hydrophobic_aa = NULL;
	int n_hydrophobic = 0;
	for (i = 1; i < chain->NAA; i++) {
	    if (chain->aa[i].id == 'A' ||
		chain->aa[i].id == 'C' || chain->aa[i].id == 'I' ||
		chain->aa[i].id == 'L' || chain->aa[i].id == 'M' ||
		chain->aa[i].id == 'F' || chain->aa[i].id == 'W' ||
		chain->aa[i].id == 'V' || chain->aa[i].id == 'p') {
		n_hydrophobic++;
		hydrophobic_aa = (int *)realloc(hydrophobic_aa,n_hydrophobic*sizeof(int));
		hydrophobic_aa[n_hydrophobic-1] = i;
	    }
	}

	w = 1.0 / (n_hydrophobic);

	for (i = 0; i < n_hydrophobic; i++)
		add(rc, rc, chain->aa[hydrophobic_aa[i]].ca);

	scale(rc, w, rc);

	for (i = 0; i < n_hydrophobic; i++) {
		subtract(ri, chain->aa[hydrophobic_aa[i]].ca, rc);
		d[0] += ri[0] * ri[0];
		d[1] += ri[1] * ri[1];
		d[2] += ri[2] * ri[2];
		c[0] += ri[1] * ri[2];
		c[1] += ri[2] * ri[0];
		c[2] += ri[0] * ri[1];
	}

	scale(d, w, d);
	scale(c, w, c);

	w = sqrt(d[0] + d[1] + d[2]);

	eigensmith(ri, d, c);	/* this messes with d */

	fprintf(sim_params->outfile,"Rgyr = %8.3f   PrinComps = %8.3f%8.3f%8.3f      N = %4d\n",
	       w, sqrt(ri[0]), sqrt(ri[1]), sqrt(ri[2]), n_hydrophobic);
	fprintf(sim_params->outfile,"Centroid      CA              %8.3f%8.3f%8.3f\n",
	       rc[0], rc[1], rc[2]);

	free(hydrophobic_aa);

}

/* Contact map, with contacts defined as CA-CA distances under 8 Angstroms. */
void cm_alpha_8(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{
	int i, j;

	//fprintf(stderr,"Printing contact map, NAA %d\n",chain->NAA);
	//if there is a contact, print 1, and 0 otherwise
	for (i = 1; i < chain->NAA; i++) {
	    for (j = 1; j < chain->NAA; j++) {
		if (distance(chain->aa[i].ca,chain->aa[j].ca)<64.0) {
		    fputc('1',sim_params->outfile);
		} else {
		    fputc('0',sim_params->outfile);
		}
	        if (j!=chain->NAA-1) fputc(' ',sim_params->outfile);
	    }
	    fputc('\n',sim_params->outfile);
	}

}

/* Print a contact map with all CB-CB contacts within the contact cutoff. */
//TODO: multi-chain protein
void cm_native_go(Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm)
{

	//allocate memory if needed
	if (biasmap->distb == NULL) {
		(biasmap)->distb = (double *) realloc((biasmap)->distb, chain->NAA * chain->NAA * sizeof(double));
		(biasmap)->NAA = chain->NAA;
		if ((biasmap)->distb == NULL) {
			stop("cm_native_go: Insufficient memory");
		}
	}

	//calculate contact map according to the actual contacts (1: contact, 0: no contact)
	fprintf(stderr,"Calculating CB--CB contact map, using a cutoff of %g.\n",sqrt(sim_params->protein_model.touch2));
	for (int i = 1; i < chain->NAA; i++) {
		Distb(i, i) = 0;
		for (int j = 1; j < i - 1; j++) {
			Distb(i, j) = 0;
			Distb(j, i) = 0;
			if (contact((chain->aa) + i, (chain->aa) + j, &(sim_params->protein_model)) && (chain->aa)[i].id != 'G' && (chain->aa)[j].id != 'G') {
				Distb(i, j) = 1;
				Distb(j, i) = 1;
			}
		}
	}
	//print contact map (X: contact, O: no contact)
	for (int i = 1; i < chain->NAA; i++) {
		for (int j = 1; j < chain->NAA; j++) {
			if (Distb(i,j)==1) {
				fputc('X',sim_params->outfile);
			} else if (Distb(i,j)==0) {
				fputc('O',sim_params->outfile);
			} else {
				stop("Unknown value in contact map.");
			}
		}
		fputc('\n',sim_params->outfile);
	}
}

/* Print fasta sequence */
void fasta(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm)
{

	for (int i = 1; i < chain->NAA; i++) {
	    fputc(chain->aa[i].id,sim_params->outfile);
	}
	fputc('\n',sim_params->outfile);

}

/* calculating Hbond pattern an individual snapshot */
void hbond_pattern (Chain *chain, Biasmap *biasmap,simulation_params *sim_params, void *mpi_comm) {

	unsigned int donor, accep;
	AA *a, *b;

	int *forward_pattern = NULL;
	int *backward_pattern = NULL;

	forward_pattern = (int *)malloc(chain->NAA * sizeof(int));
	backward_pattern = (int *)malloc(chain->NAA * sizeof(int));
	for (int i = 0; i< chain->NAA; i++) {
		forward_pattern[i] = 0;
		backward_pattern[i] = 0;
	}

	for (int i = 5; i < chain->NAA - 4; i++) {
		donor = accep = 0;
		a = chain->aa + i;

		for (int j = 1; j < chain->NAA; j++) {
			if (i == j)
				continue;

			b = (chain->aa) + j;

			/* H-pattern within the same chain only */
			if (a->chainid != b->chainid)
				continue;

			if (hdonor(a, b, &(sim_params->protein_model))) {
				//fprintf(stderr,"Found hbond don: %d to acc: %d\n",i,j);
				donor |= 1 << (15 - j + i);
				if (i<j) forward_pattern[j-i]++;
				else backward_pattern[i-j]++;
			}
			if (hdonor(b, a, &(sim_params->protein_model))) {
				//fprintf(stderr,"Found hbond don: %d to acc: %d\n",j,i);
				accep |= 1 << (15 - i + j);
				if (i<j) backward_pattern[j-i]++;
				else forward_pattern[i-j]++;
			}
		}
		//printf("%3d   %.8X   %.8X\n", i, donor, accep);
		//printf("%3d   %8X   %8X\n", i, donor, accep);
	}
	//printf("Backward_pattern ");
	for (int i = chain->NAA -1; i > 0; i--) {
		fprintf(sim_params->outfile,"%d ", backward_pattern[i]);
	}
		fprintf(sim_params->outfile,"0");
	//printf("Forward_pattern ");
	for (int i = 1; i < chain->NAA; i++) {
		fprintf(sim_params->outfile," %d", forward_pattern[i]);
	}
	fprintf(sim_params->outfile,"\n");

	free(forward_pattern);
	free(backward_pattern);

}
