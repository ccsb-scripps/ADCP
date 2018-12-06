/*
** This is a core of Metropolis Monte Carlo sampling procedure for
** simplified polypeptides. It contains the energy functions and
** the simulation procedure.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

/* MC move and metropolis criteria */
void transmutate(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params);
int flipChain(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params);
int rotate_cyclic(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params);
int transopt(Chain * chain, Chaint *chaint, Biasmap *biasmap, double ampl, double logLstar, double * currE, simulation_params *sim_params, int mod);
int move(Chain *chain, Chaint *chaint, Biasmap *biasmap,double logLstar, double *currE,int changeamp, simulation_params *sim_params);
void finalize(Chain *chain, Chaint *chaint, Biasmap *biasmap);
