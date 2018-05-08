/*
** This is a collection of functions to manipulate matrices and triplets,
** and to calculate rotation matrices and Euler angles.
**
** Copyright (c) 2004-2010 Alexei Podtelezhnikov
*/

/*
********** Data types
*/

typedef double matrix[3][3];
typedef vector triplet[3];

/*
** Triplet and matrix manipulations
*/

void casttriplet(triplet, triplet);
void transset(matrix, triplet);
void matrixvector(vector, matrix, vector);
void vectortriplet(vector, vector, triplet);
void rotation(triplet, matrix, triplet);
void fixtriplet(triplet);
void printout(triplet);
void randvector(vector);

/*
** Rotation matrix and Euler angles
*/

double *sphereframe(vector, triplet, double, double, double);
void rotmatrix(matrix, vector, double);
void eulerset(triplet, double, double, double);
double euler_bend(triplet, triplet);
double euler_twist(triplet, triplet);
double euler_alpha(triplet, triplet);
double euler_beta(triplet, triplet);
double euler_gamma(triplet, triplet);
void tripletcmp(double *, double *, triplet, triplet);

/*
** Complex number multiplication with phase accumulation (2D-rotation)
*/

struct phasor {
	double y, x;
	int k;
};

struct phasor phasiply(struct phasor, struct phasor);
int rephase(struct phasor *);
double phase(struct phasor);
