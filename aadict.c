/*
** Amino acid lexicon and conversions.
** This module contains everything about the protein models
** that can be used in the program, including all side chain
** properties.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nik Burkoff, Csilla Varnai and David Wild
*/

#include<stdio.h>
#include<ctype.h>
#include<float.h>
#include<stdlib.h>
#include<math.h>
#include"error.h"
#include"params.h"
#include"aadict.h"
#include"vector.h" /* PI/180 */


/***********************************************************/
/****  AMINO ACID CHARACTER, CODE AND INDEX CONVERSION  ****/
/***********************************************************/


/* this array contains standard 3-letter amino acid codes
the entry indicies correspond to 1-letter codes mod 32 */
const char *aac[] = {
	"XAA", "ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY",
	"HIS", "ILE", "XLE", "LYS", "LEU", "MET", "ASN", "PYL",
	"PRO", "GLN", "ARG", "SER", "THR", "SEC", "VAL", "TRP",
	"XAA", "TYR", "GLX", "XAA", "XAA", "XAA", "XAA", "XAA"
};

/* this array contains standard 1-letter amino acid codes
unused letters are 0-ed out */
const char aa_id[] = {
	'0', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
	'H', 'I', '0', 'K', 'L', 'M', 'N', 'O',
	'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
	'X', 'Y', 'Z', '0', '0', '0', '0', '0'
};

/* post-transcriptionally modified amino acids, etc */
const char *unu[] = { "MSE:M", "SAH:C", "PCA:E", "SEP:S", "TPO:T", "PTR:Y",
	"TYS:Y", "KCX:K", "LLP:K", "CSO:C", "CSD:C", "CME:C", "HYP:P", "MLY:K",
	"M3L:K", "OCS:C", "CSW:C", "CSX:C", "CSS:C", "ALY:K", "CGU:E", "CXM:M",
	"HIC:H", "SMC:C", "MLZ:K", "MEN:N", "NEP:H", "CSE:C", "SME:M", "PHD:D",
	"AGM:R", "NIY:Y", "LYZ:K", "P1L:C", "CSP:C", "MEQ:Q"
};


/* Return the amino acid one letter id, given its three letter code. */
char aa321(char *cod)
{
	int i;
	const char *fit;

	for (i = 1; i < 27; i++) {	/* looking for a letter */
		fit = aac[i];
		if (cod[0] == fit[0] && cod[1] == fit[1] && cod[2] == fit[2])
			return (i | 0x40);	/* make it a character */
	}

	for (i = 0; i < sizeof(unu) / sizeof(unu[0]); i++) {
		fit = unu[i];
		if (cod[0] == fit[0] && cod[1] == fit[1] && cod[2] == fit[2])
			return fit[4];
	}

	return 'X';
}

/* Function converts id -> index */
inline int convert_to_index(char id){
  int ans = id - 'A' + 1;
  if(ans < 1 || ans > 26) return -1;
  return ans;	
}


/***********************************************************/
/****         AMINO ACID  SIDE CHAIN PROPERTIES         ****/
/****     LIBRARY INITIALISATION AND QUERY FUNCTIONS    ****/
/***********************************************************/

/* Setting all side chain property values for an amino acid */
void initialize_one_sidechain_properties(sidechain_properties_ *sidechain_properties,
					 int num, char id, double number_of_gamma,
					 double beta_gamma1_distance,  double alpha_beta_gamma1_angle,
					 double beta_gamma2_distance,  double alpha_beta_gamma2_angle,
					 double sidechain_vdw_radius1,  double sidechain_vdw_radius2,
					 double sidechain_vdw_depth1,  double sidechain_vdw_depth2,
					 double sidechain_dihedral_gauche_plus_prob, double sidechain_dihedral_gauche_minus_prob,
					 int hydrophobic_atoms, double hydrophobic_contact_radius_CB, double hydrophobic_contact_radius_G1, double hydrophobic_contact_radius_G2,
					 int charged_atom, double charged_atom_charge,
					 int hydrogen_bond_donor_atoms, int hydrogen_bond_acceptor_atoms, double hydrogen_bond_donor_radius, double hydrogen_bond_acceptor_radius ){
  sidechain_properties[num].id = id;
  sidechain_properties[num].beta_gamma1_distance = beta_gamma1_distance;
  sidechain_properties[num].alpha_beta_gamma1_angle = alpha_beta_gamma1_angle;
  sidechain_properties[num].beta_gamma2_distance = beta_gamma2_distance;
  sidechain_properties[num].alpha_beta_gamma2_angle = alpha_beta_gamma2_angle;
  sidechain_properties[num].sidechain_vdw_radius1 = sidechain_vdw_radius1;
  sidechain_properties[num].sidechain_vdw_radius2 = sidechain_vdw_radius2;
  sidechain_properties[num].sidechain_vdw_depth1 = sidechain_vdw_depth1;
  sidechain_properties[num].sidechain_vdw_depth2 = sidechain_vdw_depth2;
  if(sidechain_vdw_depth1 < 0){
	  sidechain_properties[num].sidechain_vdw_depth1_sqrt = sidechain_vdw_depth1;
  }
  else{
    sidechain_properties[num].sidechain_vdw_depth1_sqrt = sqrt(sidechain_vdw_depth1);
  }
  if(sidechain_vdw_depth2 < 0){
	  sidechain_properties[num].sidechain_vdw_depth2_sqrt = sidechain_vdw_depth2;
  }
  else{
    sidechain_properties[num].sidechain_vdw_depth2_sqrt = sqrt(sidechain_vdw_depth2);
  }
   
  sidechain_properties[num].sidechain_dihedral_gauche_plus_prob = sidechain_dihedral_gauche_plus_prob;
  sidechain_properties[num].sidechain_dihedral_gauche_minus_prob = sidechain_dihedral_gauche_minus_prob;
  sidechain_properties[num].hydrophobic_atoms = hydrophobic_atoms;
  sidechain_properties[num].hydrophobic_contact_radius_CB = hydrophobic_contact_radius_CB;
  sidechain_properties[num].hydrophobic_contact_radius_G1 = hydrophobic_contact_radius_G1;
  sidechain_properties[num].hydrophobic_contact_radius_G2 = hydrophobic_contact_radius_G2;
  sidechain_properties[num].charged_atom = charged_atom;
  sidechain_properties[num].charged_atom_charge = charged_atom_charge;
  sidechain_properties[num].hydrogen_bond_donor_atoms = hydrogen_bond_donor_atoms;
  sidechain_properties[num].hydrogen_bond_acceptor_atoms = hydrogen_bond_acceptor_atoms;
  sidechain_properties[num].hydrogen_bond_donor_radius = hydrogen_bond_donor_radius;
  sidechain_properties[num].hydrogen_bond_acceptor_radius = hydrogen_bond_acceptor_radius;
}

/* Setting all side chain property values for all amino acids.
   This routine contains all the default paramaeters. */
void initialize_sidechain_properties(model_params *mod_params){
/*************************************************************************/
/* R. Srinivasan et al., PNAS 96(25), 14258--14263 (1999)                */
/*    + modifications:                                                   */
/*       1) added electrostatic interactions for ASP, GLU, LYS and ARG   */
/*       2) PRO geometry (CB--CG bond length and CA--CB--CG angle)       */
/*          is taken from Ho et al., Prot. Sci. 14(4) 1011--1018 (2005)  */
/*       3) side chains can take trans, gauche+ and gauche- values       */
/*          with probabilities found in the ASTRAL-1.75/abcd database    */
/*       4) added more H-bond acceptor side chains                       */
/*       5) added H-bond donor side chains                               */
/*       6) vdW interactions are different, we use LJ                    */
/*       7) PRO made amphipathic                                         */
/*       8) added hydrophobic-polar interactions between the LINUS       */
/*          hydrophobic residues' hydrophobic atoms and the polar        */
/*          residues' (B,D,E,H,K,N,Q,R,S,Z) gamma atoms                  */
/*    + potential (not yet added) modifications:                         */
/*       9) change minimum sequence separation for hydrophobic and       */
/*          electrostatic interactions from i,i+2                        */
/*      10) add partial charges or dipoles for polar side chains         */
/*************************************************************************/
/*  CORRECT_GAMMA: place gamma atoms where they actually are             */
/*  CORRECT_KMQR_GAMMA: place gamma atoms of K,M,Q and R where they      */
/*       actually are, otherwise use LINUS                               */
/*************************************************************************/
    double rs = mod_params->rs;
    double ro = mod_params->ro;
    double rcb = mod_params->rcb;
    //double rring = mod_params->rring;
	double rr = 3.0;
	double rlc = 3.5;
	double rsc = 2.5;

    double es = mod_params->vdw_depth_s;
    double eo = mod_params->vdw_depth_o;
    double ecb = mod_params->vdw_depth_cb;
    /*ro was 1.35, rs 1.8, rcb 1.65 */

    if (mod_params->use_gamma_atoms == LINUS_GAMMA) {
    /*                                                                   num  id numG d(CB,G) angle(CA,CB,G)   d(CB,G2) angle(CA,CB,G2)   RvdW_G RvdW_G2 EvdW_G EvdW_G2 Pchi1(g+) Pchi1(g-) h.phobic      Rhph_CB Rhph_G Rhph_G2  el.static q(e)   H-don H-acc    RHB_d  RHB_a */
    //fprintf(stderr,"CRANKITE GAMMA\n");
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 0,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 1,  'A', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        CB_,             2.0,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    /* Use N for B (no electrostatics)*/                                                                                    
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 2,  'B', 0,  2.130, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0,       0,        G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 3,  'C', 1,  1.822, 114.4 * M_PI_180,     -1,     -1,               rs,    -1,    es,     -1,  0.2153,  0.4890,   CB_ | G__,       2.0,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    /*Nik has made D G_ hydrophobic so charge-hydrophobic contact penalty can be used if wanted */                             
    //initialize_one_sidechain_properties(mod_params->sidechain_properties, 4,  'D', 1,  2.130, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1665,  0.5260,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 4,  'D', 1,  2.130, 112.6 * M_PI_180,     -1,     -1,              rsc,    -1,   ecb,     -1,  0.1665,  0.5260,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  2.25  );

	//initialize_one_sidechain_properties(mod_params->sidechain_properties, 5,  'E', 1,  3.000, 130.0 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0889,  0.5702,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  3.50  );
	initialize_one_sidechain_properties(mod_params->sidechain_properties, 5,  'E', 1,  3.000, 130.0 * M_PI_180,     -1,     -1,              rlc,    -1,   ecb,     -1,  0.0889,  0.5702,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  3.50  );    
	
	//initialize_one_sidechain_properties(mod_params->sidechain_properties, 6,  'F', 1,  2.908, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0983,  0.5553,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 6,  'F', 1,  2.908, 113.8 * M_PI_180,     -1,     -1,              rr,     -1,   ecb,     -1,  0.0983,  0.5553,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );

	initialize_one_sidechain_properties(mod_params->sidechain_properties, 7,  'G', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );

    //initialize_one_sidechain_properties(mod_params->sidechain_properties, 8,  'H', 1,  2.665, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1301,  0.5378,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,  2.25  );
	initialize_one_sidechain_properties(mod_params->sidechain_properties, 8,  'H', 1,  2.665, 113.8 * M_PI_180,     -1,     -1,              rr,     -1,   ecb,     -1,  0.1301,  0.5378,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,  2.25  );
    
	initialize_one_sidechain_properties(mod_params->sidechain_properties, 9,  'I', 2,  2.176, 110.5 * M_PI_180,  1.530,  110.4 * M_PI_180,   rcb,   rcb,   ecb,    ecb,  0.1402,  0.7700,   CB_ | G__ | G2_, 2.0,  2.00,  2.0,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,10,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );

    //initialize_one_sidechain_properties(mod_params->sidechain_properties,11,  'K', 1,  4.700, 130.0 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0784,  0.5582,   G__,              -1,  2.00,   -1,     G__,      1,        G__,  0x0,  4.25,    -1  );
	initialize_one_sidechain_properties(mod_params->sidechain_properties,11,  'K', 1,  4.700, 130.0 * M_PI_180,     -1,     -1,              rlc,    -1,   ecb,     -1,  0.0784,  0.5582,   G__,              -1,  2.00,   -1,     G__,      1,        G__,  0x0,  4.25,    -1  );
	
	//initialize_one_sidechain_properties(mod_params->sidechain_properties,12,  'L', 1,  2.176, 116.3 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0119,  0.6574,   CB_ | G__,       2.0,  3.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
	initialize_one_sidechain_properties(mod_params->sidechain_properties,12,  'L', 1,  2.176, 116.3 * M_PI_180,     -1,     -1,              rsc,    -1,   ecb,     -1,  0.0119,  0.6574,   CB_ | G__,       2.0,  3.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );

	
	//initialize_one_sidechain_properties(mod_params->sidechain_properties,13,  'M', 1,  3.532, 129.4 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0748,  0.6068,   CB_ | G__,       2.0,  3.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    //initialize_one_sidechain_properties(mod_params->sidechain_properties,14,  'N', 1,  2.130, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1290,  0.5605,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  3.00,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,13,  'M', 1,  3.532, 129.4 * M_PI_180,     -1,     -1,              rlc,    -1,   ecb,     -1,  0.0748,  0.6068,   CB_ | G__,       2.0,  3.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,14,  'N', 1,  2.130, 112.6 * M_PI_180,     -1,     -1,              rsc,    -1,   ecb,     -1,  0.1290,  0.5605,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  3.00,  2.25  );

    initialize_one_sidechain_properties(mod_params->sidechain_properties,15,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    /* Proline's probabilities are for +- 30 not 60 */                                                                         
    initialize_one_sidechain_properties(mod_params->sidechain_properties,16,  'P', 1,  1.501, 103.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.4602,  0.5397,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );

    //initialize_one_sidechain_properties(mod_params->sidechain_properties,17,  'Q', 1,  3.000, 130.0 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0729,  0.6052,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  4.25,  3.50  );
    //initialize_one_sidechain_properties(mod_params->sidechain_properties,18,  'R', 1,  4.900, 134.1 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0893,  0.5791,   G__,              -1,  2.00,   -1,     G__,      1,        G__,  0x0,  4.25,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,17,  'Q', 1,  3.000, 130.0 * M_PI_180,     -1,     -1,              rlc,    -1,   ecb,     -1,  0.0729,  0.6052,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  4.25,  3.50  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,18,  'R', 1,  4.900, 134.1 * M_PI_180,     -1,     -1,              rlc,    -1,   ecb,     -1,  0.0893,  0.5791,   G__,              -1,  2.00,   -1,     G__,      1,        G__,  0x0,  4.25,    -1  );
	
	initialize_one_sidechain_properties(mod_params->sidechain_properties,19,  'S', 1,  1.417, 111.1 * M_PI_180,     -1,     -1,               ro,    -1,    eo,     -1,  0.4583,  0.3010,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  2.00,  1.30  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,20,  'T', 2,  1.433, 109.9 * M_PI_180,  1.521,  109.3 * M_PI_180,    ro,   rcb,    eo,    ecb,  0.4794,  0.4418,   G2_,              -1,    -1,  2.0,     0x0,      0,        G__,  G__,  2.00,  1.30  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,21,  'U', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,22,  'V', 2,  1.521, 110.5 * M_PI_180,  1.521,  110.5 * M_PI_180,   rcb,   rcb,   ecb,    ecb,  0.0705,  0.2013,   CB_ | G__ | G2_, 2.0,  2.00,  2.0,     0x0,      0,        0x0,  0x0,    -1,    -1  );

    //initialize_one_sidechain_properties(mod_params->sidechain_properties,23,  'W', 1,  2.908, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1422,  0.4457,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,  3.00,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,23,  'W', 1,  2.908, 113.8 * M_PI_180,     -1,     -1,               rr,    -1,   ecb,     -1,  0.1422,  0.4457,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,  3.00,    -1  );

    initialize_one_sidechain_properties(mod_params->sidechain_properties,24,  'X', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );

    //initialize_one_sidechain_properties(mod_params->sidechain_properties,25,  'Y', 1,  3.314, 113.9 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1133,  0.5352,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
	initialize_one_sidechain_properties(mod_params->sidechain_properties,25,  'Y', 1,  3.314, 113.9 * M_PI_180,     -1,     -1,               rr,    -1,   ecb,     -1,  0.1133,  0.5352,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );

	
	/*Use E for Z*/                                                                                                         
    initialize_one_sidechain_properties(mod_params->sidechain_properties,26,  'Z', 1,  3.000, 130.0 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0889,  0.5702,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,27,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,28,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,29,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,30,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  ); 
    
    } else if (mod_params->use_gamma_atoms == CORRECT_KMQR_GAMMA) {
  //fprintf(stderr,"Setting up amino acid constants for CORRECT_KMQR_GAMMA model.\n");
  //crankite distances with original distances for KMQR
    /*                                                                   num  id numG d(CB,G) angle(CA,CB,G)   d(CB,G2) angle(CA,CB,G2)   RvdW_G RvdW_G2 EvdW_G EvdW_G2 Pchi1(g+) Pchi1(g-) h.phobic      Rhph_CB Rhph_G Rhph_G2  el.static q(e)   H-don H-acc    RHB_d  RHB_a */
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 0,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 1,  'A', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        CB_,             2.0,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    // Use N for B (no electrostatics)                                                                                         
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 2,  'B', 0,  2.130, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0,       0,        G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 3,  'C', 1,  1.822, 114.4 * M_PI_180,     -1,     -1,               rs,    -1,    es,     -1,  0.2153,  0.4890,   CB_ | G__,       2.0,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    // Nik has made D G_ hydrophobic so charge-hydrophobic contact penalty can be used if wanted                               
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 4,  'D', 1,  2.130, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1665,  0.5260,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 5,  'E', 1,  3.000, 130.0 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0889,  0.5702,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  3.50  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 6,  'F', 1,  2.908, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0983,  0.5553,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 7,  'G', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 8,  'H', 1,  2.665, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1301,  0.5378,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 9,  'I', 2,  2.176, 110.5 * M_PI_180,  1.530,  110.4 * M_PI_180,   rcb,   rcb,   ecb,    ecb,  0.1402,  0.7700,   CB_ | G__ | G2_, 2.0,  2.00,  2.0,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,10,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,11,  'K', 1,  1.520, 113.7 * M_PI_180,	   -1,     -1,	            rcb,    -1,   ecb,     -1,  0.0784,  0.5582,   G__,              -1,  2.00,   -1,     G__,	     1,	       G__,  0x0,  4.25,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,12,  'L', 1,  2.176, 116.3 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0119,  0.6574,   CB_ | G__,       2.0,  3.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,13,  'M', 1,  1.520, 114.5 * M_PI_180,	   -1,     -1,	            rcb,    -1,   ecb,     -1,  0.0748,  0.6068,   CB_ | G__,	    2.0,  3.00,   -1,     0x0,      0,	       0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,14,  'N', 1,  2.130, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1290,  0.5605,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  3.00,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,15,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    // Proline's probabilities are for +- 30 not 60                                                                            
    initialize_one_sidechain_properties(mod_params->sidechain_properties,16,  'P', 1,  1.501, 103.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.4602,  0.5397,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,17,  'Q', 1,  1.530, 113.5 * M_PI_180,     -1,     -1,	            rcb,    -1,   ecb,     -1,  0.0729,  0.6052,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  4.25,  3.50  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,18,  'R', 1,  1.520, 113.8 * M_PI_180,     -1,     -1,	            rcb,    -1,   ecb,     -1,  0.0893,  0.5791,   G__,              -1,  2.00,   -1,     G__,      1,        G__,  0x0,  4.25,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,19,  'S', 1,  1.417, 111.1 * M_PI_180,     -1,     -1,               ro,    -1,    eo,     -1,  0.4583,  0.3010,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  2.00,  1.30  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,20,  'T', 2,  1.433, 109.9 * M_PI_180,  1.521,  109.3 * M_PI_180,    ro,   rcb,    eo,    ecb,  0.4794,  0.4418,   G2_,              -1,    -1,  2.0,     0x0,      0,        G__,  G__,  2.00,  1.30  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,21,  'U', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,22,  'V', 2,  1.521, 110.5 * M_PI_180,  1.521,  110.5 * M_PI_180,   rcb,   rcb,   ecb,    ecb,  0.0705,  0.2013,   CB_ | G__ | G2_, 2.0,  2.00,  2.0,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,23,  'W', 1,  2.908, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1422,  0.4457,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,  3.00,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,24,  'X', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,25,  'Y', 1,  3.314, 113.9 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1133,  0.5352,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    // Use E for Z                                                                                                          
    initialize_one_sidechain_properties(mod_params->sidechain_properties,26,  'Z', 1,  3.000, 130.0 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0889,  0.5702,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,27,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,28,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,29,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,30,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  ); 

    } else if (mod_params->use_gamma_atoms == CORRECT_GAMMA) {
    /*Using actual gamma positions */
    /*                                                                   num  id numG d(CB,G) angle(CA,CB,G)   d(CB,G2) angle(CA,CB,G2)   RvdW_G RvdW_G2 EvdW_G EvdW_G2 Pchi1(g+) Pchi1(g-) h.phobic      Rhph_CB Rhph_G Rhph_G2  el.static q(e)   H-don H-acc    RHB_d  RHB_a */
    //fprintf(stderr,"PLACE GAMMA\n");
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 0,  '0', 0,     -1,    -1,                -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 1,  'A', 0,     -1,    -1,                -1,	   -1,               -1,    -1,    -1,     -1,  0,       0,        CB_,             2.0,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    /* Use N for B (no electrostatics)*/                		       				                                                                        
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 2,  'B', 0,   1.52, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0,       0,        G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 3,  'C', 1,   1.82, 113.8 * M_PI_180,     -1,     -1,               rs,    -1,    es,     -1,  0.2153,  0.4890,   CB_ | G__,       2.0,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    /*Nik has made D G_ hydrophobic so charge-hydropcontact penalty can be used if wanted 			  hobic   cont act penalty can be used if wanted */                     
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 4,  'D', 1,   1.52, 113.2 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1665,  0.5260,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 5,  'E', 1,   1.52, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0889,  0.5702,   G__,              -1,  2.00,   -1,     G__,     -1,        0x0,  G__,    -1,  3.50  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 6,  'F', 1,   1.51, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0983,  0.5553,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 7,  'G', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 8,  'H', 1,    1.5, 113.5 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1301,  0.5378,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties, 9,  'I', 2,   1.53, 110.0 * M_PI_180,   1.53,  110.0 * M_PI_180,   rcb,   rcb,   ecb,    ecb,  0.1402,  0.7700,   CB_ | G__ | G2_, 2.0,  2.00,  2.0,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,10,  '0', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,11,  'K', 1,   1.52, 113.7 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0784,  0.5582,   G__,              -1,  2.00,   -1,     G__,      1,        G__,  0x0,  4.25,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,12,  'L', 1,   1.53, 116.3 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0119,  0.6574,   CB_ | G__,       2.0,  3.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,13,  'M', 1,   1.52, 114.5 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0748,  0.6068,   CB_ | G__,       2.0,  3.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,14,  'N', 1,   1.52, 112.6 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1290,  0.5605,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  3.00,  2.25  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,15,  '0', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    /* Proline's probabilities are for +- 30 not 60      		         				  */                                                                    
    initialize_one_sidechain_properties(mod_params->sidechain_properties,16,  'P', 1,    1.5, 104.0 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.4602,  0.5397,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,17,  'Q', 1,   1.53, 113.5 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0729,  0.6052,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  4.25,  3.50  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,18,  'R', 1,   1.52, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0893,  0.5791,   G__,              -1,  2.00,   -1,     G__,      1,        G__,  0x0,  4.25,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,19,  'S', 1,   1.41, 110.8 * M_PI_180,     -1,     -1,               ro,    -1,    eo,     -1,  0.4583,  0.3010,   G__,              -1,  2.00,   -1,     0x0,      0,        G__,  G__,  2.00,  1.30  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,20,  'T', 2,   1.43, 109.1 * M_PI_180,   1.53,  111.6 * M_PI_180,    ro,   rcb,    eo,    ecb,  0.4794,  0.4418,   G2_,              -1,    -1,  2.0,     0x0,      0,        G__,  G__,  2.00,  1.30  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,21,  'U', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,22,  'V', 2,   1.53, 110.6 * M_PI_180,   1.52,  110.5 * M_PI_180,   rcb,   rcb,   ecb,    ecb,  0.0705,  0.2013,   CB_ | G__ | G2_, 2.0,  2.00,  2.0,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,23,  'W', 1,    1.5, 113.7 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1422,  0.4457,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,  3.00,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,24,  'X', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,25,  'Y', 1,   1.51, 113.7 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.1133,  0.5352,   CB_ | G__,       2.0,  3.25,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    /*Use E for Z*/                                    		        				                                                                        
    initialize_one_sidechain_properties(mod_params->sidechain_properties,26,  'Z', 1,   1.52, 113.8 * M_PI_180,     -1,     -1,              rcb,    -1,   ecb,     -1,  0.0889,  0.5702,   G__,              -1,  2.00,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,27,  '0', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,28,  '0', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,29,  '0', 0,     -1,   	-1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    initialize_one_sidechain_properties(mod_params->sidechain_properties,30,  '0', 0,     -1,    -1,	               -1,     -1,               -1,    -1,    -1,     -1,  0,       0,        0x0,              -1,    -1,   -1,     0x0,      0,        0x0,  0x0,    -1,    -1  );
    
    } else if (mod_params->use_gamma_atoms != NO_GAMMA) {
	char error_string[DEFAULT_LONG_STRING_LENGTH]="";
	sprintf(error_string,"Unknown value for use_gamma_atoms (%d)",mod_params->use_gamma_atoms);
	stop(error_string);
    }
	  
  return;
}
 


/* Function to return the CB-G distance and CA-CB-G angle for an amino acid
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   which_gamma:  1 or 2, referring to the two possible gamma atoms
   r: CB-G distance
   theta: CA-CB-G angle */
int beta_gamma_dist(char id, int which_gamma, double *r, double *theta, sidechain_properties_ *sidechain_properties) {

	int i;
	if((i = convert_to_index(id)) != -1){
		if (which_gamma == 1) {
		    *r = sidechain_properties[i].beta_gamma1_distance;
		    *theta = sidechain_properties[i].alpha_beta_gamma1_angle;
		    return 0;
		} 
		else if (which_gamma == 2) {
		    *r = sidechain_properties[i].beta_gamma2_distance;
		    *theta = sidechain_properties[i].alpha_beta_gamma2_angle;
		    return 0;
		} 
		else {
		    fprintf(stderr,"Unknown value %d of which_gamma, must be 1 or 2.\n",which_gamma);
		    return 1;
		}
    }
	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
    return 0;
}

/* Function to return the N-CA-CB-CG sidechain dihedral angle for an amino acid
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   only 1st gamma atom is needed, peptide works out where 2nd one should go */
double sidechain_dihedral(char id, sidechain_properties_ *sidechain_properties) {

	int i;
	double p_plus60 = 1/3.;
	double p_minus60 = 1/3.;
	double ans;
	if((i = convert_to_index(id)) != -1){
	  p_plus60 = sidechain_properties[i].sidechain_dihedral_gauche_plus_prob;
	  p_minus60 = sidechain_properties[i].sidechain_dihedral_gauche_minus_prob;
	  double u = rand()/(double)RAND_MAX;
	  if(u < p_plus60) ans = 60 * M_PI_180;
	  else if(u < p_plus60 + p_minus60) ans =  -60 * M_PI_180;
	  else ans = 180 * M_PI_180; 
	  //Proline is special, it has +- 30 not +- 60
	  if(id == 'P') ans /= 2;
	  return ans;
	}

	fprintf(stderr,"Error in calculating side chain dihedral angles for residue %c\n",id);
    return 1;  
}   

/* Function to return the N-CA-CB-G2 sidechain dihedral angle for an amino acid
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   only 1st gamma atom is needed, peptide works out where 2nd one should go */
double sidechain_dihedral2(char id, double chi, sidechain_properties_ *sidechain_properties){
  double chi2 = DBL_MAX;
  if(id == 'T' || id == 'I'){
    if(fabs(chi - M_PI_180 *180) < 1e-5) chi2 = 60 * M_PI_180;    	  
    if(fabs(chi - M_PI_180 * 60) < 1e-5) chi2 = -60 * M_PI_180;  
    if(fabs(chi - M_PI_180 * -60) < 1e-5) chi2 = 180 * M_PI_180;  
  }
  else if(id == 'V'){
	if(fabs(chi - M_PI_180 *180) < 1e-5) chi2 = -60 * M_PI_180;    	  
    if(fabs(chi - M_PI_180 * 60) < 1e-5) chi2 = 180 * M_PI_180;  
    if(fabs(chi - M_PI_180 * -60) < 1e-5) chi2 = 60 * M_PI_180;    
  }
  //if(chi2 == 0) fprintf(stderr,"Error: cannot calculate chi2 for id %c\n",id);
  return chi2;
}

/* Function to return the most common sidechain vdw radius for an amino acid
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   which_gamma:  1 or 2, referring to the two possible gamma atoms */
double sidechain_vdw_radius(char id, int which_gamma, sidechain_properties_ *sidechain_properties) {

	int i;

	if((i = convert_to_index(id)) != -1){
		if (which_gamma == 1) {
			 return sidechain_properties[i].sidechain_vdw_radius1;
		} 
		else if (which_gamma == 2) {
			 return sidechain_properties[i].sidechain_vdw_radius2;
		} 
		else {
			fprintf(stderr,"Unknown value %d of which_gamma, must be 1 or 2.\n",which_gamma);
		}
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}

/* Function to return the most common sidechain vdw depth for an amino acid
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   which_gamma:  1 or 2, referring to the two possible gamma atoms */
double sidechain_vdw_depth(char id, int which_gamma, sidechain_properties_ *sidechain_properties) {

	int i;

	if((i = convert_to_index(id)) != -1){
		if (which_gamma == 1) {
			 return sidechain_properties[i].sidechain_vdw_depth1;
		} 
		else if (which_gamma == 2) {
			 return sidechain_properties[i].sidechain_vdw_depth2;
		} 
		else {
			fprintf(stderr,"Unknown value %d of which_gamma, must be 1 or 2.\n",which_gamma);
		}
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}


/* Function to return the most common sidechain vdw depth sqrted for an amino acid
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   which_gamma:  1 or 2, referring to the two possible gamma atoms */
double sidechain_vdw_depth_sqrt(char id, int which_gamma, sidechain_properties_ *sidechain_properties) {

	int i;

	if((i = convert_to_index(id)) != -1){
		if (which_gamma == 1) {
			 return sidechain_properties[i].sidechain_vdw_depth1_sqrt;
		} 
		else if (which_gamma == 2) {
			 return sidechain_properties[i].sidechain_vdw_depth2_sqrt;
		} 
		else {
			fprintf(stderr,"Unknown value %d of which_gamma, must be 1 or 2.\n",which_gamma);
		}
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}

/* Function to return the hydrophobic atoms
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid */
unsigned int hydrophobic_atoms_list(char id, sidechain_properties_ *sidechain_properties) {

	int i;
	if((i = convert_to_index(id)) != -1){ 
	    return sidechain_properties[i].hydrophobic_atoms;
	}
	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}

/* Function to return the hydrophobic contact radius for an amino acid's CB, G or G2 atom
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   atom:  binary atom code */
double hydrophobic_contact_radius(char id, int atom, sidechain_properties_ *sidechain_properties) {

	int i;
	if((i = convert_to_index(id)) != -1){ 
		if (atom == CB_) {
		    return sidechain_properties[i].hydrophobic_contact_radius_CB;
		} 
		else if (atom == G__) {
		    return sidechain_properties[i].hydrophobic_contact_radius_G1;
		} 
		else if (atom == G2_) {
			return sidechain_properties[i].hydrophobic_contact_radius_G2;
		}
		else {
			fprintf(stderr,"Unknown value %x of binary atom code, must be %x, %x or %x.\n",atom, CB_, G__, G2_);
		}
	}
	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}

/* Function to return the charge of the side chain bead
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid */
double charge(char id, sidechain_properties_ *sidechain_properties) {

	int i;

	if((i = convert_to_index(id)) != -1){ 
	    return sidechain_properties[i].charged_atom_charge;
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}

/* Function to return whether a certain side chain atom is a H-bond donor
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   atom:  the binary code of the side chain atom (e.g. G__) */
int hbond_donor(char id, int atom, sidechain_properties_ *sidechain_properties) {

	int i;

	if ((i = convert_to_index(id)) != -1){
	   if  (sidechain_properties[i].hydrogen_bond_donor_atoms & atom) {
	      return 1;
	   } else {
	      return 0;
	   }
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0;
}

/* Function to return whether a certain side chain atom is a H-bond acceptor
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid
   atom:  the binary code of the side chain atom (e.g. G__) */
int hbond_acceptor(char id, int atom, sidechain_properties_ *sidechain_properties) {

	int i;

	if ((i = convert_to_index(id)) != -1){
	   if  (sidechain_properties[i].hydrogen_bond_acceptor_atoms & atom) {
	      return 1;
	   } else {
	      return 0;
	   }
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0;
}

/* Function to return the side chain radius as an H-bond donor
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid */
double sidechain_hbond_donor_radius(char id, sidechain_properties_ *sidechain_properties) {

	int i;

	if((i = convert_to_index(id)) != -1){
	   return sidechain_properties[i].hydrogen_bond_donor_radius;
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}

/* Function to return the side chain radius as an H-bond acceptor
   that was stored in the sidechain_properties library.
   id: 1-letter code of the amino acid */
double sidechain_hbond_acceptor_radius(char id, sidechain_properties_ *sidechain_properties) {

	int i;

	if((i = convert_to_index(id)) != -1){
	   return sidechain_properties[i].hydrogen_bond_acceptor_radius;
	}

	fprintf(stderr, "invalid amino acid character %c%d%d\n",id,isalpha(id),isupper(id));
	return 0.0;
}
