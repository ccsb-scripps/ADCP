/*
** This is a collection of functions to manipulate matrices and triplets,
** and to calculate rotation matrices and Euler angles.
**
** Copyright (c) 2004-2010 Alexei Podtelezhnikov
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"vector.h"
#include"rotation.h"

/* intentionally dangerous internal macro */
#ifndef WITH_SINCOS
#define sin_cos(si, co, x)	si = sin(x); co = cos(x)
#else
#define sin_cos(si, co, x)	asm ("fsincos" : "=t" (co), "=u" (si) : "0" (x))
#endif

void casttriplet(triplet x, triplet x0)
{
	x[0][0] = x0[0][0];
	x[0][1] = x0[0][1];
	x[0][2] = x0[0][2];
	x[1][0] = x0[1][0];
	x[1][1] = x0[1][1];
	x[1][2] = x0[1][2];
	x[2][0] = x0[2][0];
	x[2][1] = x0[2][1];
	x[2][2] = x0[2][2];
}

/* t=trans(x) */
void transset(matrix t, triplet x)
{
	t[0][0] = x[0][0];
	t[0][1] = x[1][0];
	t[0][2] = x[2][0];
	t[1][0] = x[0][1];
	t[1][1] = x[1][1];
	t[1][2] = x[2][1];
	t[2][0] = x[0][2];
	t[2][1] = x[1][2];
	t[2][2] = x[2][2];
}

void matrixvector(vector x, matrix t, vector y)
{
	x[0] = dotprod(t[0], y);
	x[1] = dotprod(t[1], y);
	x[2] = dotprod(t[2], y);
}

void vectortriplet(vector x, vector y, triplet t)
{
	x[0] = y[0] * t[0][0] + y[1] * t[1][0] + y[2] * t[2][0];
	x[1] = y[0] * t[0][1] + y[1] * t[1][1] + y[2] * t[2][1];
	x[2] = y[0] * t[0][2] + y[1] * t[1][2] + y[2] * t[2][2];
}

void rotation(triplet x, matrix t, triplet x0)
{
	matrixvector(x[2], t, x0[2]);
	matrixvector(x[1], t, x0[1]);
	crossprod(x[0], x[1], x[2]);
}

void fixtriplet(triplet x)
{
	double c1;

	normalize_1(x[2]);
	c1 = dotprod(x[1], x[2]);
	fling(x[1], x[1], -c1, x[2]);
	normalize_1(x[1]);
	crossprod(x[0], x[1], x[2]);
}

void printout(triplet x)
{
	printf("(%g %g %g) %g\n", x[0][0], x[0][1], x[0][2], square(x[0]));
	printf("(%g %g %g) %g\n", x[1][0], x[1][1], x[1][2], square(x[1]));
	printf("(%g %g %g) %g\n", x[2][0], x[2][1], x[2][2], square(x[2]));
	printf(" %g %g %g  %g\n",
	       dotprod(x[0], x[1]), dotprod(x[1], x[2]), dotprod(x[2], x[0]),
	       triprod(x[0], x[1], x[2]));
}

/* a vector uniformly distributed on a sphere (Knop, 1970; Marsaglia, 1972) */
void randvector(vector z)
{
	double x1, x2, ll, l;
	const double discrete = 2.0 / RAND_MAX;

	do {
		x1 = discrete * rand() - 1.0;
		x2 = discrete * rand() - 1.0;
		ll = x1 * x1 + x2 * x2;
	} while (ll > 1.0);

	z[0] = 1.0 - 2.0 * ll;

	l = 2.0 * sqrt(1.0 - ll);
	z[1] = x1 * l;
	z[2] = x2 * l;
}

/* from local spherical coordinates to lab cartesian coordinates */
double *sphereframe(vector a, triplet x, double r, double theta, double chi)
{
	double st, ct, sc, cc;
	vector b;

	sin_cos(st, ct, theta);
	sin_cos(sc, cc, chi);

	b[0] = r * st * cc;
	b[1] = r * st * sc;
	b[2] = r * ct;

	vectortriplet(a, b, x);

	return a;
}

/* setting rotation matrix given the direction and the angle
using the well-known quaternion-to-matrix transformation */
void rotmatrix(matrix t, vector x, double alpha)
{
	double si, co, q0, q1, q2, q3, p;

	alpha *= 0.5;

	sin_cos(si, co, alpha);

	si *= M_SQRT2;
	co *= M_SQRT2;

	q0 = co;
	q1 = x[0] * si;
	q2 = x[1] * si;
	q3 = x[2] * si;

	p = q0 * q0 - 1.0;
	t[0][0] = p + q1 * q1;
	t[1][1] = p + q2 * q2;
	t[2][2] = p + q3 * q3;

	t[0][1] = t[1][0] = q1 * q2;
	t[1][2] = t[2][1] = q2 * q3;
	t[2][0] = t[0][2] = q3 * q1;

	p = q0 * q1;
	t[2][1] += p;
	t[1][2] -= p;

	p = q0 * q2;
	t[0][2] += p;
	t[2][0] -= p;

	p = q0 * q3;
	t[1][0] += p;
	t[0][1] -= p;
}

/* setting triplet given the three Euler angles (floating x-convention) */
void eulerset(triplet x, double alpha, double beta, double gamma)
{
	double sa, ca, sb, cb, sc, cc;

	sin_cos(sa, ca, alpha);
	sin_cos(sb, cb, beta);
	sin_cos(sc, cc, gamma);

	x[0][0] = ca * cc - sa * cb * sc;
	x[0][1] = sa * cc + ca * cb * sc;
	x[0][2] = sb * sc;

	x[1][0] = -ca * sc - sa * cb * cc;
	x[1][1] = -sa * sc + ca * cb * cc;
	x[1][2] = sb * cc;

	x[2][0] = sa * sb;
	x[2][1] = -ca * sb;
	x[2][2] = cb;
}

/* the Euler beta angle (fail-safe, no gimbal lock) */
double euler_bend(triplet x, triplet x0)
{
	double cosb, sinbx, sinby;

	sinbx = dotprod(x[2], x0[0]);
	sinby = dotprod(x[2], x0[1]);
	cosb = dotprod(x[2], x0[2]);

	return atan2(sqrt(sinbx * sinbx + sinby * sinby), cosb);
}

/* the Euler alpha+gamma angle (fail-safe, no gimbal lock) */
double euler_twist(triplet x, triplet x0)
{
	double pcosac, psinac;

	pcosac = dotprod(x[0], x0[0]) + dotprod(x[1], x0[1]);
	psinac = dotprod(x[0], x0[1]) - dotprod(x[1], x0[0]);

	return atan2(psinac, pcosac);
}

/* the next three functions recover Euler angles away from the gimbal lock */
double euler_alpha(triplet x, triplet x0)
{
	return atan2(dotprod(x[2], x0[0]), -dotprod(x[2], x0[1]));
}

double euler_beta(triplet x, triplet x0)
{
	return acos(dotprod(x[2], x0[2]));
}

double euler_gamma(triplet x, triplet x0)
{
	return atan2(dotprod(x0[2], x[0]), dotprod(x0[2], x[1]));
}

/* fast comparison of triplets in terms of convenient measures */
void tripletcmp(double *halfsqbend, double *halfsqtwist, triplet x, triplet x0)
{
	double cosb, cosac;

	cosb = dotprod(x[2], x0[2]);
	*halfsqbend = 1.0 - cosb;

	cosac = (dotprod(x[0], x0[0]) + dotprod(x[1], x0[1])) / (1.0 + cosb);
	*halfsqtwist = 1.0 - cosac;
}

/*
** Complex number multiplication with phase accumulation (2D-rotation)
*/

struct phasor phasiply(struct phasor a, struct phasor b)
{
	struct phasor c;

	c.x = a.x * b.x - a.y * b.y;
	c.y = a.x * b.y + a.y * b.x;
	c.k = a.k + b.k;

	if (a.y > 0.0 && b.y > 0.0 && c.y < 0.0)
		c.k++;
	if (a.y < 0.0 && b.y < 0.0 && c.y > 0.0)
		c.k--;

	return c;
}

/* reduce overflowing phasor */
int rephase(struct phasor *a)
{
	int ye, xe, ee;

	ye = ilogb(a->y);
	xe = ilogb(a->x);

	ee = -((ye > xe) ? ye : xe);

	a->y = scalbn(a->y, ee);
	a->x = scalbn(a->x, ee);

	return ee;
}

double phase(struct phasor a)
{
	return atan2(a.y, a.x) + 2.0 * M_PI * a.k;
}
