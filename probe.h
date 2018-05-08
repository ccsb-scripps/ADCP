/*
** These routines calculate and print important polypeptide characteristics.
**
** Copyright (c) 2004 - 2007 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#define BEFORE_INIT 0x1;
#define AFTER_INIT 0x10;
#define BEFORE_AND_AFTER_INIT 0x11;

void tests(Chain *chain, Biasmap *biasmap, unsigned int, simulation_params *sim_params, int init_mask, void *mpi_comm);
void helps(void);

void CA_geometry(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void initialize_displacement(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void pdbout(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void sstructure(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void phipsi(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void hpattern(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void hbtot(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void hbond_geometry(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void stepsize(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void ergtot(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void n_contacts(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void n_native_contacts(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void rgyr(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void hydrophobic_rgyr(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void cm_txt(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void cm_pbm(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void cm_ideal(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void test_flex(Chain *chain, Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void cm_ideal_4(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void cm_native_go(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void hbss(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void writhe(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void torsion(Chain *chain,Biasmap *biasmap, double*,double*,int);
void energy_contributions(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void energy_gradient_wrt_parameters(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void exclude_energy_contributions(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void cos_dihedral_naac(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void strand_bias_distances(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void hydrophobic_distances(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void number_of_contacts(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void atomic_distances(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void vdw_max_gamma(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void vdw_contributions(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void checkpoint_out(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void hbond_pattern(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
//Nested Sampling only
void evidence(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void information(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void cm_alpha_8(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
void fasta(Chain *chain,Biasmap *biasmap, simulation_params *sim_params, void *mpi_comm);
