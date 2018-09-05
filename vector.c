/*
** This is a collection of procedures to determine vector sums and products,
** as well as common angles between vectors and their trigonometric functions.
**
** Copyright (c) 2003-2010 Alexei Podtelezhnikov
*/

#include<math.h>
#include"vector.h"

/*
********** Trivial vector functions
*/

double *castvec(vector b, vector a)
{
	b[0] = a[0];
	b[1] = a[1];
	b[2] = a[2];

	return b;
}

/* vector addition */
double *add(vector apb, vector a, vector b)
{
	apb[0] = a[0] + b[0];
	apb[1] = a[1] + b[1];
	apb[2] = a[2] + b[2];

	return apb;
}

/* vector subtraction */
double *subtract(vector amb, vector a, vector b)
{
	amb[0] = a[0] - b[0];
	amb[1] = a[1] - b[1];
	amb[2] = a[2] - b[2];

	return amb;
}

/* scale a vector */
double *scale(vector qa, scalar q, vector a)
{
	qa[0] = q * a[0];
	qa[1] = q * a[1];
	qa[2] = q * a[2];

	return qa;
}

/* update a vector */
double *fling(vector apb, vector a, double bb, vector b)
{
	apb[0] = a[0] + bb * b[0];
	apb[1] = a[1] + bb * b[1];
	apb[2] = a[2] + bb * b[2];

	return apb;
}

/* linear combination */
double *lincomb(vector apb, double aa, vector a, double bb, vector b)
{
	apb[0] = aa * a[0] + bb * b[0];
	apb[1] = aa * a[1] + bb * b[1];
	apb[2] = aa * a[2] + bb * b[2];

	return apb;
}

/* dot product */
double dotprod(vector a, vector b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* square length of a vector */
double square(vector a)
{
	return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

/* invsquare calculates inverse square length 1/(aa) */
double invsquare(vector a)
{
	return 1.0 / square(a);
}

/* cross product */
double *crossprod(vector ab, vector a, vector b)
{
	ab[0] = a[1] * b[2] - a[2] * b[1];
	ab[1] = a[2] * b[0] - a[0] * b[2];
	ab[2] = a[0] * b[1] - a[1] * b[0];

	return ab;
}

/* entry-wise Schur-Hadamard product */
double *schurprod(vector ab, vector a, vector b)
{
	ab[0] = a[0] * b[0];
	ab[1] = a[1] * b[1];
	ab[2] = a[2] * b[2];

	return ab;
}

/* scalar triple product */
double triprod(vector a, vector b, vector c)
{
	vector bc;

	crossprod(bc, b, c);

	return dotprod(a, bc);
}

/* vector triangle area */
double *triarea(vector abc, vector a, vector b, vector c)
{
	vector ab, bc;

	subtract(ab, a, b);
	subtract(bc, b, c);
	crossprod(abc, ab, bc);

	return abc;
}

/* square distance between two points */
double distance(vector a, vector b)
{
	double x, y, z;

	x = a[0] - b[0];
	y = a[1] - b[1];
	z = a[2] - b[2];

	return x * x + y * y + z * z;
}

/* square point-line distance */
double pointline(vector a, vector b, vector v)
{
	vector ab, abv;

	subtract(ab, a, b);
	crossprod(abv, ab, v);

	return square(abv) / square(v);
}

/* square line-line distance */
double lineline(vector a, vector b, vector u, vector v)
{
	vector ab, uv;
	double vol;

	subtract(ab, a, b);
	crossprod(uv, u, v);
	vol = dotprod(ab, uv);

	return vol * vol / square(uv);
}

/* non-orthogonal projected components */
double *triprjct(vector q, vector p, vector a, vector b, vector c)
{
	vector abc, x;
	double d;

	triarea(abc, a, b, c);
	d = 1.0 / square(abc);

	triarea(x, p, b, c);
	q[0] = dotprod(x, abc) * d;
	triarea(x, a, p, c);
	q[1] = dotprod(x, abc) * d;
	triarea(x, a, b, p);
	q[2] = dotprod(x, abc) * d;

	return q;
}

/* non-orthogonal 2d components a = q0*b + q1*c */
double *twocomps(vector q, vector a, vector b, vector c)
{
	double ab, bc, ca, b2, c2, d;

	ab = dotprod(a, b);
	bc = dotprod(b, c);
	ca = dotprod(c, a);

	b2 = square(b);
	c2 = square(c);

	d = 1.0 / (b2 * c2 - bc * bc);

	q[0] = (ab * c2 - bc * ca) * d;
	q[1] = (b2 * ca - ab * bc) * d;
	q[2] = 1.0 - q[0] - q[1];

	return q;
}

/* non-orthogonal 3d components */
double *tricomps(vector q, vector p, vector a, vector b, vector c)
{
	double d;
	vector x;

	crossprod(x, b, c);
	d = 1.0 / dotprod(a, x);

	q[0] = dotprod(p, x) * d;
	crossprod(x, c, a);
	q[1] = dotprod(p, x) * d;
	crossprod(x, a, b);
	q[2] = dotprod(p, x) * d;

	return q;
}

/* vector normalization */
double normalize(vector a)
{
	double inva;

	inva = 1.0 / sqrt(square(a));
	scale(a, inva, a);

	return inva;
}

/* quick normalization of an almost unit vector */
double normalize_1(vector a)
{
	double inva, a2;

	a2 = square(a);
	/* Taylor expansion */
	inva = 1.5 - 0.5 * a2;
	/* inva = 1.875 - (1.25 - 0.375 * a2) * a2; */
	/* inva = 2.1875 - (2.1875 - (1.3125 - 0.3125 * a2) * a2) * a2; */
	/* inva = 0.5 * (1.0 / a2 + 1.0); */
	/* Newton-Halley iteration */
	/* inva = (a2 + 3.0) / (3.0 * a2 + 1.0); */

	scale(a, inva, a);

	return inva;
}

/*
********** Angles between vectors (Arcfunctions)
*/

/* xy_ functions facilitate implicit summation of angles */
Comp xy_add(Comp a, Comp b)
{
	Comp xy;

	xy.x = a.x * b.x - a.y * b.y;
	xy.y = a.x * b.y + a.y * b.x;

	return xy;
}

Comp xy_angle(vector a, vector b)
{
	vector ab;
	Comp xy;

	xy.x = dotprod(a, b);
	crossprod(ab, a, b);
	xy.y = sqrt(square(ab));

	return xy;
}

Comp xy_dihedral(vector a, vector b, vector c)
{
	vector ab, bc;
	Comp xy;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	xy.x = dotprod(ab, bc);
	xy.y = dotprod(ab, c) * sqrt(square(b));

	return xy;
}

/* The angle between two vectors is given by
	angle = atan { |[ab]| / (ab) }
Computing atan2 is significantly faster than acos. */
double angle(vector a, vector b)
{
	vector ab;

	crossprod(ab, a, b);

	return atan2(sqrt(square(ab)), dotprod(a, b));
}

/* dihedral calculates the angle given by three vectors
	dihedral = atan { [abc]|b| / ([ab][bc]) }
Using atan2 saves CPU cycles and defines the proper quadrant.
Total: 22 multiplications, 12 additions, 1 sqrt and 1 atan2. */
double dihedral(vector a, vector b, vector c)
{
	vector ab, bc;
	double b1;

	crossprod(ab, a, b);
	crossprod(bc, b, c);
	b1 = sqrt(square(b));

	return atan2(dotprod(ab, c) * b1, dotprod(ab, bc));
}

/* If |b| = 1, dihedral angle calculation avoids sqrt. It is very fast. */
double dihedral_1(vector a, vector b, vector c)
{
	vector ab, bc;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	return atan2(dotprod(ab, c), dotprod(ab, bc));
}

double dihedral_4(vector a0, vector a1, vector a2, vector a3)
{
	vector a, b, c;

	subtract(a, a1, a0);
	subtract(b, a2, a1);
	subtract(c, a3, a2);

	return dihedral(a, b, c);
}

double dihedral_rama(vector a0, vector a1, vector a2, vector a3, double b1)
{
	vector a, b, c;

	subtract(a, a1, a0);
	subtract(b, a2, a1);
	subtract(c, a3, a2);

	vector ab, bc;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	return atan2(dotprod(ab, c) * b1, dotprod(ab, bc));
}

/* solid angle between three vectors according to Oosterom and Strackee (1983).
Total: 33 multiplications, 20 additions, 3 sqrt's, and 1 atan2  */
double excess(vector a, vector b, vector c)
{
	double abc, ab, bc, ca, a1, b1, c1;

	abc = triprod(a, b, c);
	ab = dotprod(a, b);
	bc = dotprod(b, c);
	ca = dotprod(c, a);
	a1 = sqrt(square(a));
	b1 = sqrt(square(b));
	c1 = sqrt(square(c));

	return 2.0 * atan2(abc, a1 * b1 * c1 + ab * c1 + bc * a1 + ca * b1);
}

/*
********** Direct trigonometric functions from vectors
*/

/* Square cosine of the angle between two vectors.
This is the fastest angular measure */
double sqcosine(vector a, vector b)
{
	double ab, a2, b2;

	ab = dotprod(a, b);
	a2 = square(a);
	b2 = square(b);

	return (ab * ab) / (a2 * b2);
}

/* Cosine of the angle between two vectors.
Total: 10 multiplications, 6 additions, 1 division, and 1 sqrt. */
double cosine(vector a, vector b)
{
	double ab, a2, b2;

	ab = dotprod(a, b);
	a2 = square(a);
	b2 = square(b);

	return ab / sqrt(a2 * b2);
}

/* cosine of a-b-c angle */
double cosangle(vector a, vector b, vector c) {
	vector ba;
	subtract(ba, a, b);
	vector bc;
	subtract(bc, c, b);
	return cosine(ba,bc);
}

int cosgreater(vector a, vector b, double c)
{
	double ab, a2, b2;

	ab = dotprod(a, b);
	a2 = square(a);
	b2 = square(b);

	return ab * fabs(ab) > c * fabs(c) * a2 * b2;
}

/* Sine of the angle between two vectors */
double sine(vector a, vector b)
{
	double ab, c2;
	vector c;

	ab = dotprod(a, b);
	crossprod(c, a, b);
	c2 = square(c);

	return sqrt(c2 / (ab * ab + c2));

	/* faster and less accurately
	   return sqrt(1.0 - sqcosine(a, b)); */
}

/* Sine of an angle with a phase shift specified by its y and x coordinates.
This calculates linear combination  of y * cos + x * sin. */
double phasine(vector a, vector b, double y, double x)
{
	double ab, c2;
	vector c;

	ab = dotprod(a, b);
	crossprod(c, a, b);
	c2 = square(c);

	return (y * ab + x * sqrt(c2)) / sqrt(ab * ab + c2);
}

/* Tangent of the angle between two vectors */
double tangent(vector a, vector b)
{
	vector ab;

	crossprod(ab, a, b);

	return sqrt(square(ab)) / dotprod(a, b);
}

/* Cosine of the triple angle between two vectors */
double costri(vector a, vector b)
{
	double cosab;

	cosab = cosine(a, b);

	return (4.0 * cosab * cosab - 3.0) * cosab;
}

/* Square cosine of the dihedral angle between three vectors. */
double sqcosdihedral(vector a, vector b, vector c)
{
	vector ab, bc;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	return sqcosine(ab, bc);
}

/* Cosine of dihedral angle. */
double cosdihedral(vector a, vector b, vector c)
{
	vector ab, bc;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	return cosine(ab, bc);
}

/* Sine of dihedral angle. */
double sindihedral(vector a, vector b, vector c)
{
	vector ab, bc;
	double pcos, abc, b2;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	pcos = dotprod(ab, bc);
	abc = dotprod(ab, c);
	b2 = square(b);

	return abc * sqrt(b2 / (pcos * pcos + abc * abc * b2));
}

/* Sine of dihedral angle with phase shift specified by y and x coordinates.
This calculates linear combination  of y * cos + x * sin. */
double phasindihedral(vector a, vector b, vector c, double y, double x)
{
	vector ab, bc;
	double pcos, psin;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	pcos = dotprod(ab, bc);
	psin = dotprod(ab, c) * sqrt(square(b));

	return (y * pcos + x * psin) / sqrt(pcos * pcos + psin * psin);
}

/* Tangent of dihedral angle */
double tandihedral(vector a, vector b, vector c)
{
	vector ab, bc;
	double pcos, psin;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	pcos = dotprod(ab, bc);
	psin = dotprod(ab, c) * sqrt(square(b));

	return psin / pcos;
}

/* Cosine of triple dihedral angle */
double costridihedral(vector a, vector b, vector c)
{
	vector ab, bc;

	crossprod(ab, a, b);
	crossprod(bc, b, c);

	return costri(ab, bc);
}
