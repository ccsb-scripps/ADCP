/*
 * flex.c
 *
 *  Created on: 14 Sep 2012
 *      Author: nik
 */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<float.h>
#include<math.h>
#include<sys/stat.h>

#include"error.h"
#include"params.h"
#include"vector.h"
#include"rotation.h"
#include"aadict.h"
#include"peptide.h"
#include"vdw.h"
#include"flex.h"
#include"energy.h"
#include"probe.h"
#include"checkpoint_io.h"
#include"metropolis.h"

#ifdef PARALLEL
#include<mpi.h>
#endif

double ca_ha_distance = 1.085;
double cb_hb_distance = 1.085;
double ca_cb_hb_angle = 109.5*M_PI_180;

double sg_hg_distance = 1.32;
double og_hg_distance = 0.95;

/*To check below */
double cb_sg_hg_angle = 108*M_PI_180;
double cb_og_hg_angle = 108*M_PI_180;
double ca_cb_og_hg_dihedral = 180; //do not need M_PI_180
double ca_cb_sg_hg_dihedral = 180; //do not need M_PI_180


int output_one_aa_covin(AA *aa,int N_atom,int *C_atom, int use_gamma_atoms, FILE*fptr){
  int bond_number = 5;
  fprintf(fptr,"%d %d %d\n",N_atom,N_atom+1,bond_number+(aa->id=='P'? 1 : 0)); //N -CA
  fprintf(fptr,"%d %d %d\n",N_atom+1,N_atom+2,bond_number); //CA-C
  *C_atom = N_atom+2;
  fprintf(fptr,"%d %d %d\n",N_atom+2,N_atom+3,bond_number); //C-O
  if(aa->id != 'P'){
    fprintf(fptr,"%d %d %d\n",N_atom,N_atom+4,bond_number); //N-H
    fprintf(fptr,"%d %d %d\n",N_atom+1,N_atom+5,bond_number); //CA - CB(or H for Glycine)
    fprintf(fptr,"%d %d %d\n",N_atom+1,N_atom+6,bond_number); //CA - H
    if(aa->id == 'G') return N_atom+7;
    fprintf(fptr,"%d %d %d\n",N_atom+5,N_atom+7,bond_number); //CB - C/O/SG (or H for alanine)
    fprintf(fptr,"%d %d %d\n",N_atom+5,N_atom+8,bond_number); //CB - CG2 (or H)
    fprintf(fptr,"%d %d %d\n",N_atom+5,N_atom+9,bond_number); //CB - H
    if(aa->id == 'A' || use_gamma_atoms == 0 ) return N_atom+10;
    switch (aa->id) {
    case 'I':
    case 'V':
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+10,bond_number); //CG1-HG11
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+11,bond_number); //CG1-HG12
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+12,bond_number); //CG1-HG13
      fprintf(fptr,"%d %d %d\n",N_atom+8,N_atom+13,bond_number); //CG2-HG21
      fprintf(fptr,"%d %d %d\n",N_atom+8,N_atom+14,bond_number); //CG2-HG22
      fprintf(fptr,"%d %d %d\n",N_atom+8,N_atom+15,bond_number); //CG2-HG23
      return N_atom+16;
    case 'T':
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+10,bond_number); //OG-HG1
      fprintf(fptr,"%d %d %d\n",N_atom+8,N_atom+11,bond_number); //CG2-HG21
      fprintf(fptr,"%d %d %d\n",N_atom+8,N_atom+12,bond_number); //CG2-HG22
      fprintf(fptr,"%d %d %d\n",N_atom+8,N_atom+13,bond_number); //CG2-HG23
      return N_atom+14;
    case 'S':
    case 'C':
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+10,bond_number); // O/S - H
      return N_atom+11;

    default:
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+10,bond_number); //CG-HG1
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+11,bond_number); //CG-HG2
      fprintf(fptr,"%d %d %d\n",N_atom+7,N_atom+12,bond_number); //CG-HG3
      return N_atom+13;

    }//end switch
  }
  //proline here
  fprintf(fptr,"%d %d %d\n",N_atom+1,N_atom+4,bond_number); //CA - CB
  fprintf(fptr,"%d %d %d\n",N_atom+1,N_atom+5,bond_number); //CA - H
  fprintf(fptr,"%d %d %d\n",N_atom+4,N_atom+6,bond_number); //CB - CG
  fprintf(fptr,"%d %d %d\n",N_atom+4,N_atom+7,bond_number); //CB - H2
  fprintf(fptr,"%d %d %d\n",N_atom+4,N_atom+8,bond_number); //CB - H3
  if(use_gamma_atoms == 0) return N_atom+9;
#ifdef LINUS_1995
  return N_atom+9; //No CG on LINUS PRO
#endif
  fprintf(fptr,"%d %d %d\n",N_atom+6,N_atom+9,bond_number); //CG-HG1
  fprintf(fptr,"%d %d %d\n",N_atom+6,N_atom+10,bond_number); //CG-HG2
  fprintf(fptr,"%d %d %d\n",N_atom+6,N_atom+11,bond_number); //CG-HG3
  /*TODO make N_atom+9 = CD and add N_atom+11,N_atom+12 to HD atoms */
  return N_atom+12;
}


void output_covin(Chain *chain, model_params *protein_model,char *filename){
  FILE *fptr;
  if((fptr = fopen(filename, "r"))){
    fclose(fptr);
    fprintf(stderr,"Note overwriting %s\n",filename);
  }
  fptr = fopen(filename,"w");

  int N_atom = 1;
  int C_atom = 0;
  int i;

  chain->flex_data->oxy_index[0] = -1;

  for(i = 1; i < chain->NAA; i++){
    chain->flex_data->oxy_index[i] = N_atom + 3;
    N_atom = output_one_aa_covin(&(chain->aa[i]),N_atom,&C_atom,protein_model->use_gamma_atoms,fptr);
    if(i != chain->NAA -1){
      fprintf(fptr,"%d %d 6\n",C_atom,N_atom); //peptide_bond
    }
  }

  fclose(fptr);

}

void output_stackedin(char* filename){
  FILE *fptr;
  if((fptr = fopen(filename, "r"))){
    fclose(fptr);
    fprintf(stderr,"Note overwriting %s\n",filename);
  }
  fptr = fopen(filename,"w");
  fclose(fptr);
}

static void get_H_position(vector zz,vector N,vector CA, vector CB, double chi,double theta,double bondlength){
  double c1;
  chi *= M_PI_180;
  triplet x;
  vector b;

  /* local frame */
  subtract(x[2], CB, CA);
  normalize(x[2]);
  subtract(x[0], N, CA);
  c1 = dotprod(x[0], x[2]);
  fling(x[0], x[0], -c1, x[2]);
  normalize(x[0]);
  crossprod(x[1], x[2], x[0]);

  /* spherical coordinates */
  /* in the Euler coordintes PI-(ca-cb-g)angle is the theta angle */
  sphereframe(b, x, bondlength, M_PI-theta, chi);
  add(zz,CB,b);

}


static void get_HA_position(vector zz, vector CA, vector N, vector C, vector CB,int levo){
  double p,q;
  vector xy;

  vector ca_c,ca_n,ca_cb;
  subtract(ca_n,N,CA);
  subtract(ca_c,C,CA);
  subtract(ca_cb,CB,CA);
  /* normalize bond vectors */
  scale(ca_n, ican, ca_n);
  scale(ca_c, icac, ca_c);

  add(zz, ca_n, ca_c);
  crossprod(xy, ca_n, ca_c);
  q = dotprod(ca_n, ca_c);

  /* -CA- angles are flexible synchronous scissors */
  p = -.54835422911833237471;
  q = .56844147803928796138 - .75614562856111715263 * q;

  /* levo-dextro isomerism alters one sign */
  if (!levo)
      q = -q;
  lincomb(zz, p * cacb, zz, q * cacb, xy);
  normalize(zz);
  scale(zz,ca_ha_distance,zz);
  add(zz,zz,CA);
}


static int Hyd_pdbrecord( AA *a, int j, model_params *mod_params, FILE *outfile, FILE *ca_outfile){
  char fmt[] = "ATOM  %5d %4s XAA A9999    %8.3f%8.3f%8.3f\n";
  const char chain = 'A';
  char *gatom = "testingtesting";
  char *g2atom = "testingtesting";



  sprintf(fmt + 14, "%3s %c%4d", aa123(a->id), chain, a->num & 0xFFF);
  fmt[23] = ' ';

  fprintf(outfile,fmt, ++j, " N  ", a->n[0], a->n[1], a->n[2]);
  fprintf(outfile,fmt, ++j, " CA ", a->ca[0], a->ca[1], a->ca[2]);
  fprintf(ca_outfile,fmt,j," CA ",a->ca[0], a->ca[1], a->ca[2]);

  fprintf(outfile,fmt, ++j, " C  ", a->c[0], a->c[1], a->c[2]);
  fprintf(outfile,fmt, ++j, " O  ", a->o[0], a->o[1], a->o[2]);

  if(a->id == 'G'){
    fprintf(outfile,fmt, ++j, " H  ", a->h[0], a->h[1], a->h[2]);
    vector HA;
    get_HA_position(HA,a->ca,a->n,a->c,a->cb,LEV);
    fprintf(outfile,fmt, ++j, " HA2", HA[0], HA[1], HA[2]);
    get_HA_position(HA,a->ca,a->n,a->c,a->cb,!LEV);
    fprintf(outfile,fmt, ++j, " HA3", HA[0], HA[1], HA[2]);
    return j;
  }

  vector HA;
  get_HA_position(HA,a->ca,a->n,a->c,a->cb,!LEV);

  vector HB;

  if(a->id == 'P'){
    fprintf(outfile,fmt, ++j, " CB ", a->cb[0], a->cb[1], a->cb[2]);
    fprintf(outfile,fmt, ++j, " HA ", HA[0], HA[1], HA[2]);
#ifdef LINUS_1995
    get_H_position(HB,a->n,a->ca, a->cb,60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB1", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->n,a->ca, a->cb,-60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB2", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->n,a->ca, a->cb,180,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB3", HB[0], HB[1], HB[2]);
    return j;
#endif
    if(mod_params->use_gamma_atoms == NO_GAMMA){
      get_H_position(HB,a->n,a->ca, a->cb,60,ca_cb_hb_angle,cb_hb_distance);
      fprintf(outfile,fmt, ++j, " HB1", HB[0], HB[1], HB[2]);
      get_H_position(HB,a->n,a->ca, a->cb,-60,ca_cb_hb_angle,cb_hb_distance);
      fprintf(outfile,fmt, ++j, " HB2", HB[0], HB[1], HB[2]);
      get_H_position(HB,a->n,a->ca, a->cb,180,ca_cb_hb_angle,cb_hb_distance);
      fprintf(outfile,fmt, ++j, " HB3", HB[0], HB[1], HB[2]);
      return j;
    }
    gatom = " CG ";
    fprintf(outfile,fmt, ++j, gatom, a->g[0], a->g[1], a->g[2]);
    double chi[2];
    if(a->chi1 > 0){chi[0] = 145; chi[1]=-95;}
    else{chi[0] = -145; chi[1] = 95;}
    get_H_position(HB,a->n,a->ca, a->cb,chi[0],ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB2", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->n,a->ca, a->cb,chi[1],ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB3", HB[0], HB[1], HB[2]);
    /*Note have tetrahedral CH3 at CG position needs to be changed if CD included */
    get_H_position(HB,a->ca,a->cb, a->g,60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HG1", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->ca,a->cb, a->g,-60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HG2", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->ca,a->cb, a->g,180,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HG3", HB[0], HB[1], HB[2]);
    return j;
  }

  fprintf(outfile,fmt, ++j, " H  ", a->h[0], a->h[1], a->h[2]);
  fprintf(outfile,fmt, ++j, " CB ", a->cb[0], a->cb[1], a->cb[2]);
  fprintf(outfile,fmt, ++j, " HA ", HA[0], HA[1], HA[2]);
  if(mod_params->use_gamma_atoms == NO_GAMMA || a->id=='A' ){
    get_H_position(HB,a->n,a->ca, a->cb,60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB1", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->n,a->ca, a->cb,-60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB2", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->n,a->ca, a->cb,180,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB3", HB[0], HB[1], HB[2]);
    return j;
  }
  switch (a->id) {
  case 'I':
  case 'V':
    gatom = " CG1";
    g2atom = " CG2";
    break;
  case 'T':
    gatom = " OG1";
    g2atom = " CG2";
    break;
  case 'S':
    gatom = " OG ";
    break;
  case 'C':
    gatom = " SG ";
    break;
  case 'U':
    gatom = "SEG ";
    break;
  default:
    gatom = " CG ";
    break;
  }

  if(a->etc & G2_){
    fprintf(outfile,fmt, ++j, gatom, a->g[0], a->g[1], a->g[2]);
    fprintf(outfile,fmt, ++j, g2atom, a->g2[0], a->g2[1], a->g2[2]);
    double chi;
    if((a->chi1 > 120/M_180_PI || a->chi1 < 0 )&&  (a->chi2 > 120/M_180_PI || a->chi2 < 0)  ) chi = 60;
    else if((a->chi1 < -120/M_180_PI || a->chi1 > 0 )&&  (a->chi2 < -120/M_180_PI || a->chi2 > 0) ) chi = -60;
    else chi = 180;

    get_H_position(HB,a->n,a->ca, a->cb,chi,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB3", HB[0], HB[1], HB[2]);
  }
  else{
    fprintf(outfile,fmt, ++j, gatom, a->g[0], a->g[1], a->g[2]);
    double chi[2];
    if(a->chi1 < 120/M_180_PI && a->chi1 > 0){chi[0] = -60; chi[1] = 180;}
    else if(a->chi1 > -120/M_180_PI && a->chi1 < 0){chi[0] = 60; chi[1] = 180;}
    else{chi[0] = 60; chi[1] = -60;}
    get_H_position(HB,a->n,a->ca, a->cb,chi[0],ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB2", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->n,a->ca, a->cb,chi[1],ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HB3", HB[0], HB[1], HB[2]);
  }


  if(!(a->etc & G2_)){
    if(a->id=='S'){
      get_H_position(HB,a->ca,a->cb, a->g,ca_cb_og_hg_dihedral,cb_og_hg_angle,og_hg_distance);
      fprintf(outfile,fmt, ++j, " HG ", HB[0], HB[1], HB[2]);
      return j;
    }
    if(a->id =='C'){
      get_H_position(HB,a->ca,a->cb, a->g,ca_cb_sg_hg_dihedral,cb_sg_hg_angle,sg_hg_distance);
      fprintf(outfile,fmt, ++j, " HG ", HB[0], HB[1], HB[2]);
      return j;
    }
    get_H_position(HB,a->ca,a->cb, a->g,60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HG1", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->ca,a->cb, a->g,-60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HG2", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->ca,a->cb, a->g,180,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, " HG3", HB[0], HB[1], HB[2]);
    return j;
  }
  if(a->id == 'T'){
    get_H_position(HB,a->ca,a->cb, a->g,ca_cb_og_hg_dihedral,cb_og_hg_angle,og_hg_distance);
    fprintf(outfile,fmt, ++j, " OH ",HB[0], HB[1], HB[2]);
    get_H_position(HB,a->ca,a->cb, a->g2,60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, "HG21",HB[0], HB[1], HB[2]);
    get_H_position(HB,a->ca,a->cb, a->g2,-60,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, "HG22", HB[0], HB[1], HB[2]);
    get_H_position(HB,a->ca,a->cb, a->g2,180,ca_cb_hb_angle,cb_hb_distance);
    fprintf(outfile,fmt, ++j, "HG23", HB[0], HB[1], HB[2]);
    return j;
  }
  get_H_position(HB,a->ca,a->cb, a->g,60,ca_cb_hb_angle,cb_hb_distance);
  fprintf(outfile,fmt, ++j, "HG11", HB[0], HB[1], HB[2]);
  get_H_position(HB,a->ca,a->cb, a->g,-60,ca_cb_hb_angle,cb_hb_distance);
  fprintf(outfile,fmt, ++j, "HG12", HB[0], HB[1], HB[2]);
  get_H_position(HB,a->ca,a->cb, a->g,180,ca_cb_hb_angle,cb_hb_distance);
  fprintf(outfile,fmt, ++j, "HG13", HB[0], HB[1], HB[2]);
  get_H_position(HB,a->ca,a->cb, a->g2,60,ca_cb_hb_angle,cb_hb_distance);
  fprintf(outfile,fmt, ++j, "HG21", HB[0], HB[1], HB[2]);
  get_H_position(HB,a->ca,a->cb, a->g2,-60,ca_cb_hb_angle,cb_hb_distance);
  fprintf(outfile,fmt, ++j, "HG22", HB[0], HB[1], HB[2]);
  get_H_position(HB,a->ca,a->cb, a->g2,180,ca_cb_hb_angle,cb_hb_distance);
  fprintf(outfile,fmt, ++j, "HG23", HB[0], HB[1], HB[2]);

  return j;




}



/* Printing an amino acid chain in PDB format with extra hydrogens */
static void Hyd_pdbprint( AA *a, int count, model_params *mod_params, FILE *outfile, FILE *ca_outfile){
  int i, j;
  int model = 1;
  fprintf(outfile,"MODEL     %4d\n", model);

  for (j = 0, i = 1; i < count; i++)
    j = Hyd_pdbrecord(a + i, j, mod_params, outfile,ca_outfile);

  fprintf(outfile,"TER   %5d      %3s %c%4d\n",
      ++j, aa123(a[i - 1].id), 'A', i - 1);
  fprintf(outfile,"ENDMDL\n");
}


void Hyd_pdbout(Chain *chain,model_params *protein_model,char *filename, char *ca_filename){
  FILE *fptr, *fptr2;
  if((fptr = fopen(filename, "r"))){
    fclose(fptr);
    fprintf(stderr,"Note overwriting %s\n",filename);
  }
  fptr = fopen(filename,"w");
  if((fptr2 = fopen(ca_filename, "r"))){
    fclose(fptr2);
    fprintf(stderr,"Note overwriting %s\n",ca_filename);
  }
  fptr2 = fopen(ca_filename,"w");
  Hyd_pdbprint(chain->aa, chain->NAA, protein_model, fptr,fptr2);
  fclose(fptr); fclose(fptr2);
}


void output_hphobesin(Chain * chain, model_params *protein_model,char *filename, FILE *logfile){
  const double hphobe_distance_squared = 16.0;
  FILE *fptr;
  if((fptr = fopen(filename, "r"))){
    fclose(fptr);
    fprintf(stderr,"Note overwriting %s\n",filename);
  }
  fptr = fopen(filename,"w");

  int number_hphobes = 0;

  if(protein_model->kauzmann_param != 0){
    int i,j;
    for(i = 1; i < chain->NAA; i++){
      if(chain->aa[i].etc & HYDROPHOBIC){
        int which_atoms_i[3] = {0,0,0};
        if ((chain->aa[i].etc & CB_) &&  hydrophobic_contact_radius(chain->aa[i].id, CB_, protein_model->sidechain_properties) > 0 ){
          which_atoms_i[0] = chain->flex_data->oxy_index[i]+2;
        }
        if ((chain->aa[i].etc & G__) &&  hydrophobic_contact_radius(chain->aa[i].id, G__, protein_model->sidechain_properties) > 0 ){
          which_atoms_i[1] = chain->flex_data->oxy_index[i]+4;
        }
        if ((chain->aa[i].etc & G2_) &&  hydrophobic_contact_radius(chain->aa[i].id, G2_, protein_model->sidechain_properties) > 0 ){
          which_atoms_i[2] = chain->flex_data->oxy_index[i]+5;
        }

        for(j = i+protein_model->hydrophobic_min_separation; j < chain->NAA; j++){
          if(chain->aa[j].etc & HYDROPHOBIC){
            if ((chain->aa[j].etc & CB_) &&  hydrophobic_contact_radius(chain->aa[j].id, CB_, protein_model->sidechain_properties) > 0 ){
              if(which_atoms_i[0] != 0 && distance(chain->aa[i].cb,chain->aa[j].cb) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[0],chain->flex_data->oxy_index[j]+2);}
              if(which_atoms_i[1] != 0 && distance(chain->aa[i].g,chain->aa[j].cb) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[1],chain->flex_data->oxy_index[j]+2);}
              if(which_atoms_i[2] != 0 && distance(chain->aa[i].g2,chain->aa[j].cb) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[2],chain->flex_data->oxy_index[j]+2);}

            }
            if ((chain->aa[j].etc & G__) &&  hydrophobic_contact_radius(chain->aa[j].id, G__, protein_model->sidechain_properties) > 0 ){
              if(which_atoms_i[0] != 0 && distance(chain->aa[i].cb,chain->aa[j].g) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[0],chain->flex_data->oxy_index[j]+4);}
              if(which_atoms_i[1] != 0 && distance(chain->aa[i].g,chain->aa[j].g) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[1],chain->flex_data->oxy_index[j]+4);}
              if(which_atoms_i[2] != 0 && distance(chain->aa[i].g2,chain->aa[j].g) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[2],chain->flex_data->oxy_index[j]+4);}

            }
            if ((chain->aa[j].etc & G2_) &&  hydrophobic_contact_radius(chain->aa[j].id, G2_, protein_model->sidechain_properties) > 0 ){
              if(which_atoms_i[0] != 0 && distance(chain->aa[i].cb,chain->aa[j].g2) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[0],chain->flex_data->oxy_index[j]+5);}
              if(which_atoms_i[1] != 0 && distance(chain->aa[i].g,chain->aa[j].g2) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[1],chain->flex_data->oxy_index[j]+5);}
              if(which_atoms_i[2] != 0 && distance(chain->aa[i].g2,chain->aa[j].g2) < hphobe_distance_squared){number_hphobes++; fprintf(fptr,"%d %d 5\n",which_atoms_i[2],chain->flex_data->oxy_index[j]+5);}

            }
          }
        }
      }

    }


  }

  fprintf(logfile,"Number Hphobes %d ",number_hphobes);

  fclose(fptr);
}

void output_hbondsin(Chain * chain, model_params *protein_model, int use_bias_only, char *filename,double nma_hstrength, FILE *logfile){
  FILE *fptr;
  if((fptr = fopen(filename, "r"))){
    fclose(fptr);
    fprintf(stderr,"Note overwriting %s\n",filename);
  }
  fptr = fopen(filename,"w");
  int i,j;



  int number_hbonds = 0;

  if(use_bias_only ==0.0){
    if(nma_hstrength == 1.0){

      for(i = 1; i < chain->NAA; i++){
        for(j = i+2; j <chain->NAA; j++){
          if(hdonor(&(chain->aa[i]),&(chain->aa[j]),protein_model)){
            number_hbonds++;
            fprintf(fptr,"%d %d %.3lf 5\n",chain->flex_data->oxy_index[i]+1,chain->flex_data->oxy_index[j],-4.0);
          }
          if(hdonor(&(chain->aa[j]),&(chain->aa[i]),protein_model)){
            number_hbonds++;
            fprintf(fptr,"%d %d %.3lf 5\n",chain->flex_data->oxy_index[j]+1,chain->flex_data->oxy_index[i],-4.0);
          }
        }
      }
    }
    else{
      for(i = 1; i < chain->NAA; i++){
        for(j = i+2; j <chain->NAA; j++){
          if( chain->aa[i].id != 'P' && hstrength(chain->aa[i].n,  chain->aa[i].h,  chain->aa[j].o,  chain->aa[j].c, protein_model) > nma_hstrength ){
            fprintf(fptr,"%d %d %.3lf 5\n",chain->flex_data->oxy_index[i]+1,chain->flex_data->oxy_index[j],-4.0);
            number_hbonds++;
          }
          if(chain->aa[j].id != 'P' && hstrength(chain->aa[j].n,  chain->aa[j].h,  chain->aa[i].o,  chain->aa[i].c, protein_model) > nma_hstrength){
            fprintf(fptr,"%d %d %.3lf 5\n",chain->flex_data->oxy_index[j]+1,chain->flex_data->oxy_index[i],-4.0);
            number_hbonds++;
          }
        }
      }


    }
  }
  else{
    if(nma_hstrength == 1.0){

      for(i = 0; i < chain->flex_data->number_hbond; i++){
        if(hdonor(&(chain->aa[chain->flex_data->Hbond_aaH[i]]),&(chain->aa[chain->flex_data->Hbond_aaO[i]]),protein_model)){
          number_hbonds++;
          fprintf(fptr,"%d %d %.3lf 5\n",chain->flex_data->oxy_index[chain->flex_data->Hbond_aaH[i]]+1,chain->flex_data->oxy_index[chain->flex_data->Hbond_aaO[i]],-4.0);
        }
      }
    }
    else{
      for(i = 0; i < chain->flex_data->number_hbond; i++){

        if(hstrength(chain->aa[chain->flex_data->Hbond_aaH[i]].n,  chain->aa[chain->flex_data->Hbond_aaH[i]].h,  chain->aa[chain->flex_data->Hbond_aaO[i]].o,  chain->aa[chain->flex_data->Hbond_aaO[i]].c, protein_model) > nma_hstrength ){
          fprintf(fptr,"%d %d %.3lf 5\n",chain->flex_data->oxy_index[chain->flex_data->Hbond_aaH[i]]+1,chain->flex_data->oxy_index[chain->flex_data->Hbond_aaO[i]],-4.0);
          number_hbonds++;
        }


      }

    }
  }


  fprintf(logfile,"H-bonds %d ",number_hbonds);

  fclose(fptr);
}

int output_and_run_flex(Chain *chain,Biasmap *biasmap, Chain *input_chains, simulation_params *sim_params,double *rmsd){

  //if(isnan(totenergy(chain))){fprintf(stderr,"ERROR READ IN OUTPUT_AND_RUN");}

  char outpdb[DEFAULT_SHORT_STRING_LENGTH];
  strcpy(outpdb,sim_params->flex_params.output_path);
  strcat(outpdb,sim_params->flex_params.outputpdb_filename);

  char outpdbmat[DEFAULT_SHORT_STRING_LENGTH];
  strcpy(outpdbmat,sim_params->flex_params.output_path);
  strcat(outpdbmat,"pdbmat.structure");

  Hyd_pdbout(chain,&(sim_params->protein_model),outpdb,outpdbmat);


  char hbondsin[DEFAULT_SHORT_STRING_LENGTH];
  strcpy(hbondsin,sim_params->flex_params.output_path);
  strcat(hbondsin,"hbonds.in");
  output_hbondsin(chain,&(sim_params->protein_model),sim_params->flex_params.only_bias_hbonds,  hbondsin,sim_params->flex_params.hstrength_cutoff,sim_params->outfile);

  char hphobesin[DEFAULT_SHORT_STRING_LENGTH];
  strcpy(hphobesin,sim_params->flex_params.output_path);
  strcat(hphobesin,"hphobes.in");

  output_hphobesin(chain,&(sim_params->protein_model),hphobesin,sim_params->outfile);

  int k = system(sim_params->flex_params.flex_cmd);
  if (k != 0) stop("system returned non-0 exit status.");

  return read_in_after_flex(chain,biasmap,input_chains,sim_params,rmsd);
}

int read_in_after_flex(Chain *chain,Biasmap *biasmap,Chain *input_chains,simulation_params *sim_params,double *rmsd){

  int ans = 0;
  Chain *ichain = (Chain *)malloc(sizeof(Chain)); ichain->NAA = 0;
  ichain->NAA =0; ichain->aa = NULL; ichain->xaa = NULL; ichain->erg = NULL; ichain->xaa_prev = NULL;
  Chaint *ichaint= (Chaint*)malloc(sizeof(Chaint));
  ichaint->aat = NULL; ichaint->ergt = NULL;  ichaint->xaat = NULL; ichaint->xaat_prev = NULL;

  FILE *fptr2;
  char rmsd_filename[DEFAULT_LONG_STRING_LENGTH];
  sprintf(rmsd_filename,"%srmsd", sim_params->flex_params.output_path);
  fptr2 = fopen(rmsd_filename,"r");

  double *rmsd_temp = (double*)malloc(sizeof(double)*sim_params->flex_params.size_of_filename_to_read_in);
  int inchain;
  for(inchain = 0; inchain < sim_params->flex_params.size_of_filename_to_read_in; inchain++){

    FILE *fptr;


    if((fptr = fopen(sim_params->flex_params.filenames_to_read_in[inchain], "r"))){

      int k = fscanf(fptr2,"%lf;\n",&(rmsd_temp[inchain]));
      if (k != 1) stop("could not read rmsd_temp");
      chain->flex_data->read_in_flex++;

      int i; for (i = 1; pdbin(ichain,sim_params,fptr) != EOF; i++);
      fixpeptide(ichain->aa, ichain->NAA, &(sim_params->protein_model));
      chkpeptide(ichain->aa, ichain->NAA, &(sim_params->protein_model));
      update_sim_params_from_chain(ichain,sim_params); // updating NAA and seq

      initialize(ichain,ichaint,sim_params);
      energy_matrix_calculate(ichain,biasmap,&(sim_params->protein_model));

      ichain->ll = -totenergy(ichain);


      copybetween(&(input_chains[inchain]),ichain);
      ans++;
      fclose(fptr);
    }
    else{
      fscanf(fptr2,"\n");
      input_chains[inchain].ll = -DBL_MAX;

      //fprintf(stdout,"NaN ");

    }
  }

  //fprintf(stdout,"\n");

  freemem_chain(ichain);
  free(ichain);
  freemem_chaint(ichaint);
  free(ichaint);

  fclose(fptr2);

  //fprintf(stderr,"Read in %d pdbs\n",ans);
  /*Re-ordering rmsd so it is correct */
  int mode;
  int total_per_mode = sim_params->flex_params.size_of_filename_to_read_in;// (sim_params->flex_params.tonode - sim_params->flex_params.fromnode+1);
  for(mode = 0; mode <= 0 /*sim_params->flex_params.tonode-sim_params->flex_params.fromnode*/; mode++){
    int i; for(i = 0; i < total_per_mode; i++ ){
      if(i % 2  == 0) rmsd[mode*total_per_mode+i] = rmsd_temp[mode*total_per_mode+ (i/2)];
      else rmsd[mode*total_per_mode+i] = rmsd_temp[mode*total_per_mode+(total_per_mode/2)+(i-1)/2];
    }
  }
  free(rmsd_temp);

  return ans;
}

void setup_bias_hbonds(Chain *chain, Biasmap * biasmap){
  int i,j;
  /*First Beta sheets */
  for(i = 1; i < biasmap->NAA; i++){
    if(biasmap->distb[i*biasmap->NAA + i] != -1 || chain->aa[i].id == 'P') continue;
    for(j = 1; j < biasmap->NAA; j++ ){
      if(abs(j-i)<2) continue;
      if(biasmap->distb[j*biasmap->NAA + j] == -1 && biasmap->distb[i*biasmap->NAA + j] == 1.0){

          chain->flex_data->number_hbond++;
          chain->flex_data->Hbond_aaH = realloc(chain->flex_data->Hbond_aaH,sizeof(int)*chain->flex_data->number_hbond);
          chain->flex_data->Hbond_aaO = realloc(chain->flex_data->Hbond_aaO,sizeof(int)*chain->flex_data->number_hbond);
          chain->flex_data->Hbond_aaH[chain->flex_data->number_hbond-1] = i;
          chain->flex_data->Hbond_aaO[chain->flex_data->number_hbond-1] = j;
          if(j > 1){
            chain->flex_data->number_hbond++;
            chain->flex_data->Hbond_aaH = realloc(chain->flex_data->Hbond_aaH,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaO = realloc(chain->flex_data->Hbond_aaO,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaH[chain->flex_data->number_hbond-1] = i;
            chain->flex_data->Hbond_aaO[chain->flex_data->number_hbond-1] = j-1;
          }
          if(j < biasmap->NAA - 1){
            chain->flex_data->number_hbond++;
            chain->flex_data->Hbond_aaH = realloc(chain->flex_data->Hbond_aaH,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaO = realloc(chain->flex_data->Hbond_aaO,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaH[chain->flex_data->number_hbond-1] = i;
            chain->flex_data->Hbond_aaO[chain->flex_data->number_hbond-1] = j+1;
          }
          if(i>1 && chain->aa[i-1].id != 'P' && biasmap->distb[(i-1)*biasmap->NAA + (i-1) ] != -1.0){
            chain->flex_data->number_hbond++;
            chain->flex_data->Hbond_aaH = realloc(chain->flex_data->Hbond_aaH,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaO = realloc(chain->flex_data->Hbond_aaO,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaH[chain->flex_data->number_hbond-1] = i-1;
            chain->flex_data->Hbond_aaO[chain->flex_data->number_hbond-1] = j;
          }
          if(i<biasmap->NAA -1 && chain->aa[i+1].id != 'P' && biasmap->distb[(i+1)*biasmap->NAA + (i+1) ] != -1.0){
            chain->flex_data->number_hbond++;
            chain->flex_data->Hbond_aaH = realloc(chain->flex_data->Hbond_aaH,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaO = realloc(chain->flex_data->Hbond_aaO,sizeof(int)*chain->flex_data->number_hbond);
            chain->flex_data->Hbond_aaH[chain->flex_data->number_hbond-1] = i+1;
            chain->flex_data->Hbond_aaO[chain->flex_data->number_hbond-1] = j;
          }


      }

    }
  }
  /*Then alpha helix   (O)i->(H)i+4 */
  for(i = 1; i < biasmap->NAA-4; i++){
    if(biasmap->distb[i*biasmap->NAA + i]== 1.0){
      int ishelix = 1;
      for(j= i+1; j <= i+4; j++){
        if(biasmap->distb[j*biasmap->NAA + j]!= 1.0) ishelix = 0;
      }
      if(ishelix==1 && chain->aa[i+4].id != 'P'){
        chain->flex_data->number_hbond++;
        chain->flex_data->Hbond_aaH = realloc(chain->flex_data->Hbond_aaH,sizeof(int)*chain->flex_data->number_hbond);
        chain->flex_data->Hbond_aaO = realloc(chain->flex_data->Hbond_aaO,sizeof(int)*chain->flex_data->number_hbond);
        chain->flex_data->Hbond_aaH[chain->flex_data->number_hbond-1] = i+4;
        chain->flex_data->Hbond_aaO[chain->flex_data->number_hbond-1] = i;
      }
    }
  }

  //for(i = 0; i < chain->nma_data->number_hbond; i++){
  //  fprintf(stderr,"%d %d\n",chain->nma_data->Hbond_aaO[i],chain->nma_data->Hbond_aaH[i]);
  //}fflush(stderr);


}


void initialize_flex(Chain *chain, Chain **input_chains,Biasmap *biasmap, simulation_params* sim_params,double **rmsd){
  chain->flex_data = NULL;
  chain->flex_data = (FLEX_data*)malloc(sizeof(FLEX_data));

  chain->flex_data->oxy_index = NULL;
  if(!(chain->flex_data->oxy_index = malloc(sizeof(int)*chain->NAA)) ){
    stop("Error allocating memory in initialize_nma\n");
  }

  chain->flex_data->total_flex = 0;
  chain->flex_data->accepted_flex = 0;
  chain->flex_data->read_in_flex = 0;

  char create_directory[DEFAULT_SHORT_STRING_LENGTH] ;
  sprintf(create_directory,"mkdir -p %s",sim_params->flex_params.output_path);
  int k = system(create_directory);
  if (k != 0) stop("system returned non-0 exit status.");

  char covin[DEFAULT_SHORT_STRING_LENGTH];
  strcpy(covin,sim_params->flex_params.output_path);
  strcat(covin,"cov.in");
  output_covin(chain,&(sim_params->protein_model),covin);

  char stackedin[DEFAULT_SHORT_STRING_LENGTH];
  strcpy(stackedin,sim_params->flex_params.output_path);
  strcat(stackedin,"stacked.in");
  output_stackedin(stackedin);

  chain->flex_data->Hbond_aaH = NULL;
  chain->flex_data->Hbond_aaO = NULL;
  chain->flex_data->number_hbond = 0;

  if(sim_params->flex_params.only_bias_hbonds == 1){
    setup_bias_hbonds(chain,biasmap);
  }


  *input_chains=NULL;
  *input_chains = (Chain*)malloc(sizeof(Chain)*sim_params->flex_params.size_of_filename_to_read_in);
  int i;
  *rmsd = NULL;
  *rmsd = (double*)malloc(sizeof(double)*sim_params->flex_params.size_of_filename_to_read_in);

  for(i = 0; i < sim_params->flex_params.size_of_filename_to_read_in; i++){

      (*input_chains)[i].aa = NULL;
      (*input_chains)[i].NAA = chain->NAA;
      (*input_chains)[i].xaa = NULL;
      (*input_chains)[i].xaa_prev = NULL;
      (*input_chains)[i].erg = NULL;


      allocmem_chain(&((*input_chains)[i]),chain->NAA,chain->Nchains);

      int j;
      for(j = 1; j < chain->NAA; j++){
        (*input_chains)[i].aa[j].id = chain->aa[j].id;
        (*input_chains)[i].aa[j].num = chain->aa[j].num;
        (*input_chains)[i].aa[j].etc = chain->aa[j].etc;
      }

    }


}

void finalize_flex(Chain *chain,Chain **input_chains, simulation_params *sim_params, double **rmsd){
  //int i; for(i = 1; i < chain->NAA; i++) printf("%d\n",chain->nma_data.Oxy_index[i]);
  chain->flex_data->total_flex = 0;
  chain->flex_data->accepted_flex = 0;
  chain->flex_data->read_in_flex = 0;
  free(chain->flex_data->oxy_index);
  if(chain->flex_data->number_hbond !=  0){
    free(chain->flex_data->Hbond_aaH);
    free(chain->flex_data->Hbond_aaO);
    chain->flex_data->number_hbond = 0;
  }

  chain->flex_data->Hbond_aaH=NULL;
  chain->flex_data->Hbond_aaO=NULL;


  free(chain->flex_data);
  chain->flex_data = NULL;




  int i;
  for(i = 0; i < sim_params->flex_params.size_of_filename_to_read_in; i++){
    freemem_chain(&((*input_chains)[i]));
  }
  free(*input_chains);
  *input_chains = NULL;

  free(*rmsd);
  *rmsd = NULL;



}

#ifdef PARALLEL
int check_finished(simulation_params *sim_params){
  char filename[DEFAULT_SHORT_STRING_LENGTH];
  //strcpy(filename, sim_params->flex_params.output_path);
  strcpy(filename,"finished");
  struct stat buf;
  int i = stat ( filename, &buf );
  /* File found */
  if ( i == 0 ){
    char cmd[DEFAULT_SHORT_STRING_LENGTH];
    strcpy(cmd,"rm -rf ");
    strcat(cmd,filename);
    int k = system(cmd);
    if (k != 0) stop("system returned non-0 exit status.");
    strcpy(cmd,"rm -rf ");
    //strcat(cmd,sim_params->flex_params.output_path );
    strcat(cmd,"ready");
    k = system(cmd);
    if (k != 0) stop("system returned non-0 exit status.");
    return 0; //finished
  }
  return 1; //not finished

}

void update_flex_cmd_and_amplitude(int flex_iteration, simulation_params * sim_params, FLEX_data * flex_data){
  if(flex_iteration % 50 == 0){
  double current_accept_rate = flex_data->accepted_flex/ (double)flex_data->total_flex;


  if(flex_data->read_in_flex / (double)flex_data->total_flex > sim_params->flex_params.acceptance_rate_aim ){

  if(current_accept_rate > sim_params->flex_params.acceptance_rate_aim + sim_params->flex_params.acceptance_rate_tolerance && sim_params->flex_params.totalconf < 10000 ){
    sim_params->flex_params.totalconf /= sim_params->flex_params.flex_changing_factor;
    sim_params->flex_params.freq = sim_params->flex_params.totalconf;
    flex_setup_command(&(sim_params->flex_params));
    flex_param_print((sim_params->flex_params),sim_params->outfile);
  }
  else if(current_accept_rate < sim_params->flex_params.acceptance_rate_aim - sim_params->flex_params.acceptance_rate_tolerance&& sim_params->flex_params.totalconf > 75){
    sim_params->flex_params.totalconf *= sim_params->flex_params.flex_changing_factor;
    sim_params->flex_params.freq = sim_params->flex_params.totalconf;
    flex_setup_command(&(sim_params->flex_params));
    flex_param_print((sim_params->flex_params),sim_params->outfile);
  }
  }
  flex_data->read_in_flex = 0;
    flex_data->accepted_flex = 0;
    flex_data->total_flex = 0;
  }
  else{
	flex_setup_command(&(sim_params->flex_params));
  }



}

void ns_for_flex_processor(MPI_Comm FLEX_WORLD, int rank, Biasmap *biasmap, simulation_params* sim_params){
  int finished = 1;
  MPI_Status status;
  int i;
  Chain *chain = (Chain *)malloc(sizeof(Chain)); chain->NAA = 0;
  chain->NAA =0; chain->aa = NULL; chain->xaa = NULL; chain->erg = NULL; chain->xaa_prev = NULL;
  Chaint *chaint= (Chaint*)malloc(sizeof(Chaint));
  chaint->aat = NULL; chaint->ergt = NULL;  chaint->xaat = NULL; chaint->xaat_prev = NULL;

  sim_params->outfile_name = (char*)malloc(sizeof(char)*DEFAULT_SHORT_STRING_LENGTH);
  strcpy(sim_params->outfile_name,"flex.log");
  sim_params->outfile = fopen(sim_params->outfile_name,"w");
  flex_param_print(sim_params->flex_params,sim_params->outfile);

  freopen("flex.error", "w", stderr);


  MPI_Bcast(&chain->NAA,1,MPI_INT,0,FLEX_WORLD);
  allocmem_chain(chain,chain->NAA,chain->Nchains);

  for(i = 1; i < chain->NAA; i++){
    MPI_Bcast(&chain->aa[i].id,1,MPI_CHAR,0,FLEX_WORLD);
    MPI_Bcast(&chain->aa[i].num,1,MPI_INT,0,FLEX_WORLD);
    MPI_Bcast(&chain->aa[i].etc,1,MPI_INT,0,FLEX_WORLD);
  }


  update_sim_params_from_chain(chain,sim_params); // updating NAA and seq
  biasmap_initialise(chain,biasmap,&(sim_params->protein_model));
  aat_init(chain, chaint);


  Chain *input_chains;
  double *rmsd;
  initialize_flex(chain,&input_chains,biasmap,sim_params,&rmsd);


  int flex_iteration = 0;

  while(finished == 1){
    mpi_rec_chain(chain,0,rank,&sim_params->logLstar,1,FLEX_WORLD);




    flex_iteration++;
    update_flex_cmd_and_amplitude(flex_iteration,sim_params,(chain->flex_data));

    fprintf(sim_params->outfile,"FLEX Iteration %d ",flex_iteration);


    double start_ll = chain->ll;

    int number_new_chain = output_and_run_flex(chain,biasmap,input_chains,sim_params,rmsd);

    char cmd[DEFAULT_LONG_STRING_LENGTH];
    strcpy(cmd,"touch ");
    //strcat(cmd,sim_params->flex_params.output_path);
    strcat(cmd,"ready");
    int k = system(cmd);
    if (k != 0) stop("system returned non-0 exit status.");

    sprintf(cmd,"%sremove_flex.sh %s",sim_params->flex_params.flex_dir,sim_params->flex_params.output_path);
    k = system(cmd);
    if (k != 0) stop("system returned non-0 exit status.");

    finished = check_finished(sim_params);
    if(finished == 1){
      MPI_Recv(&(sim_params->logLstar),1,MPI_DOUBLE,0,10,FLEX_WORLD,&status);
      fprintf(sim_params->outfile,"Start Energy: %lf  E* %lf\n",-start_ll,-sim_params->logLstar);




      int number_to_send = 0;
      for(i = 0; i < sim_params->flex_params.size_of_filename_to_read_in; i++){
        if(input_chains[i].ll != -DBL_MAX ){
          double currE = -input_chains[i].ll;
          for(int j = 0; j <  sim_params->flex_params.MCiter /*&& currE > -logLstar*/; j++){
            move(&(input_chains[i]),chaint,biasmap,sim_params->logLstar,&currE,0, sim_params);
          }
          input_chains[i].ll = -currE;
          if(input_chains[i].ll > sim_params->logLstar) number_to_send++;
        }
      }



      fprintf(sim_params->outfile,"FLEX read_in %d  FLEX Accepted %d\n",number_new_chain,number_to_send);


      MPI_Send(&number_to_send,1,MPI_INT,0,11,FLEX_WORLD);
      int sending;
      for(sending = 0; sending < sim_params->flex_params.size_of_filename_to_read_in; sending++){
        //int which_freq = 1+(sending % (sim_params->flex_params.size_of_filename_to_read_in ))/ /*(sim_params->flex_params.tonode - sim_params->flex_params.fromnode+1)))/ */ 2;
        fprintf(sim_params->outfile,"Mode %d, %c Freq %d ",sim_params->flex_params.mode /*+ sending / (sim_params->flex_params.size_of_filename_to_read_in /(sim_params->flex_params.tonode - sim_params->flex_params.fromnode+1))*/ ,
            (sending % 2 == 0 ? '+' : '-'),/*which_freq */ sim_params->flex_params.freq );
        if(input_chains[sending].ll != -DBL_MAX){
            fprintf(sim_params->outfile," E %lf Ediff %lf RMSD %lf ", -input_chains[sending].ll,-input_chains[sending].ll+start_ll ,rmsd[sending] );
        }
        if(input_chains[sending].ll > sim_params->logLstar){
          fprintf(sim_params->outfile,"@");
            mpi_send_chain(&(input_chains[sending]), rank, 0, &sim_params->logLstar, 0, FLEX_WORLD);
        }
        fprintf(sim_params->outfile,"\n");
      }


      chain->flex_data->accepted_flex += number_to_send;
      chain->flex_data->total_flex += sim_params->flex_params.size_of_filename_to_read_in;
      fprintf(sim_params->outfile,"Overall: FLEX Accepted: %d  FLEX Attempted %d FLEX Acceptance Rate %lf \n",chain->flex_data->accepted_flex,chain->flex_data->total_flex,chain->flex_data->accepted_flex/(double)chain->flex_data->total_flex);

      fprintf(sim_params->outfile,"********************\n");
    }

  }
  freemem_chaint(chaint);
  finalize_flex(chain,&input_chains,sim_params,&rmsd);


}
#endif



int check_flex_ready(simulation_params *sim_params){
  char filename[DEFAULT_SHORT_STRING_LENGTH];
  //strcpy(filename, sim_params->flex_params.output_path);
  strcpy(filename,"ready");
  struct stat buf;
  int i = stat ( filename, &buf );

  /* File found */
  if ( i == 0 ){
    char cmd[DEFAULT_SHORT_STRING_LENGTH];
    strcpy(cmd,"rm -rf ");
    strcat(cmd,filename);
    int k = system(cmd);
    if (k != 0) stop("system returned non-0 exit status.");
    return 1; //ready
  }
  return 0; //not ready
}

