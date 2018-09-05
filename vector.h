/*
** This is a collection of procedures to determine vector sums and products,
** as well as common angles between vectors and their trigonometric functions.
**
** Copyright (c) 2003-2010 Alexei Podtelezhnikov
*/

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef M_1_PI
# define M_1_PI		0.31830988618379067154	/* 1/pi */
#endif

#ifndef M_180_PI
# define M_180_PI	57.2957795130823208768	/* 180/pi */
#endif

#ifndef M_PI_180
# define M_PI_180	.0174532925199432957692	/* pi/180 */
#endif

#ifndef M_PI_3
# define M_PI_3		1.04719755119659774615	/* pi/3 */
#endif

#ifndef M_PI_6
# define M_PI_6		0.523598775598298873075	/* pi/3 */
#endif

#ifndef M_SQRT2
# define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
#endif

#ifndef M_SQRT3
# define M_SQRT3	1.73205080756887729353	/* sqrt(3) */
#endif

/*
********** Data types
*/

typedef double scalar;
typedef double vector[3];
typedef struct {
	double x, y;
} Comp;

/*
********** Trivial vector functions
*/

double *castvec(vector, vector);
double *add(vector, vector, vector);
double *subtract(vector, vector, vector);
double *scale(vector, scalar, vector);
double *fling(vector, vector, scalar, vector);
double *lincomb(vector, scalar, vector, scalar, vector);
double dotprod(vector, vector);
double square(vector);
double invsquare(vector);
double *crossprod(vector, vector, vector);
double *schurprod(vector, vector, vector);
double triprod(vector, vector, vector);
double *triarea(vector, vector, vector, vector);
double distance(vector, vector);
double pointline(vector, vector, vector);
double lineline(vector, vector, vector, vector);
double *triprjct(vector, vector, vector, vector, vector);
double *twocomps(vector, vector, vector, vector);
double *tricomps(vector, vector, vector, vector, vector);
double normalize(vector);
double normalize_1(vector);

/*
********** Angles between vectors (Arcfunctions)
*/

Comp xy_add(Comp, Comp);
Comp xy_angle(vector, vector);
Comp xy_dihedral(vector, vector, vector);

double angle(vector, vector);
double dihedral(vector, vector, vector);
double dihedral_1(vector, vector, vector);
double dihedral_4(vector, vector, vector, vector);
double dihedral_rama(vector, vector, vector, vector, double);
double excess(vector, vector, vector);

/*
********** Trigonometric functions directly from vectors
*/

double sqcosine(vector, vector);
double cosine(vector, vector);
double cosangle(vector a,vector b, vector c);
double sine(vector, vector);
double phasine(vector, vector, double, double);
double tangent(vector, vector);
double costri(vector, vector);

double sqcosdihedral(vector, vector, vector);
double cosdihedral(vector, vector, vector);
double sindihedral(vector, vector, vector);
double phasindihedral(vector, vector, vector, double, double);
double tandihedral(vector, vector, vector);
double costridihedral(vector, vector, vector);
