/*
** This will average the input structures provided as models or separate files
** and estimate B-factor (8*pi*pi*variance).
**
** Copyright (c) 2007-2010 Alexei Podtelezhnikov
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"params.h"
#include"aadict.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"

#define VER "befa 0.1, Copyright (c) 2005-2010 Alexei Podtelezhnikov\n"
#define USE "Usage: %s [[-f] infile] [-o outfile]\n"

#define comp(a,b) a[0],a[1],a[2],1.0,bfactor(a,b)

AA *aa = NULL, *ava, *sqa;
int NAA = 0, ii;
int Nchains = 0;

double bfactor(vector x, vector x2)
{
	double bfactor;
	vector y;

	schurprod(y, x, x);
	subtract(y, x2, y);

	bfactor = 8 * M_PI * M_PI * (y[0] + y[1] + y[2]);
	if (bfactor > 999.99)
		bfactor = 999.99;

	return bfactor;
}

int pdbrec( AA *a,  AA *b, int i, int j)
{
	const char fmt[] =
	    "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
	const char *res;
	const char chain = 'A';

	res = aa123(a->id);

	printf(fmt, ++j, " N  ", res, chain, i, comp(a->n, b->n));
	printf(fmt, ++j, " CA ", res, chain, i, comp(a->ca, b->ca));
	printf(fmt, ++j, " C  ", res, chain, i, comp(a->c, b->c));
	printf(fmt, ++j, " O  ", res, chain, i, comp(a->o, b->o));
	if (a->id != 'G')
		printf(fmt, ++j, " CB ", res, chain, i, comp(a->cb, b->cb));
	if (a->id != 'P')
		printf(fmt, ++j, " H  ", res, chain, i, comp(a->h, b->h));

	return j;
}

void update(vector av, vector sq, vector x)
{
	double ri;
	vector y;

	ri = 1.0 / ii;

	subtract(y, x, av);
	scale(y, ri, y);
	add(av, av, y);

	schurprod(y, x, x);
	subtract(y, y, sq);
	scale(y, ri, y);
	add(sq, sq, y);
}

void parse_input(void)
{
	int i;

	while (getpdb(&aa, &NAA, &Nchains, stdin)) {
		ii++;
		for (i = 1; i < NAA; i++) {
			update(ava[i].h, sqa[i].h, aa[i].h);
			update(ava[i].n, sqa[i].n, aa[i].n);
			update(ava[i].ca, sqa[i].ca, aa[i].ca);
			update(ava[i].c, sqa[i].c, aa[i].c);
			update(ava[i].o, sqa[i].o, aa[i].o);
			update(ava[i].cb, sqa[i].cb, aa[i].cb);
		}
	}
}

void read_options(int argc, char *argv[])
{
	int i, opt;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			freopen(argv[i], "r", stdin);
			continue;
		}

		opt = argv[i][1];
		if (++i >= argc)
			opt = 0;

		switch (opt) {
		case 'f':
			freopen(argv[i], "r", stdin);
			break;
		case 'o':
			freopen(argv[i], "w", stdout);
			break;
		default:
			fprintf(stderr, VER USE, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char *argv[])
{
	int i, j;

	read_options(argc, argv);

	getpdb(&aa, &NAA, &Nchains, stdin);
	ii = 1;

	ava = ( AA *) malloc(sizeof(AA) * NAA);
	sqa = ( AA *) malloc(sizeof( AA) * NAA);

	memcpy(ava, aa, sizeof( AA) * NAA);

	for (i = 1; i < NAA; i++) {
		schurprod(sqa[i].h, aa[i].h, aa[i].h);
		schurprod(sqa[i].n, aa[i].n, aa[i].n);
		schurprod(sqa[i].ca, aa[i].ca, aa[i].ca);
		schurprod(sqa[i].c, aa[i].c, aa[i].c);
		schurprod(sqa[i].o, aa[i].o, aa[i].o);
		schurprod(sqa[i].cb, aa[i].cb, aa[i].cb);
		sqa[i].id = aa[i].id;
		sqa[i].num = aa[i].num;
	}

	parse_input();

	for (j = 0, i = 1; i < NAA; i++)
		j = pdbrec(ava + i, sqa + i, i, j);

	printf("TER   %5d      %3s %c%4d\n",
	       ++j, aa123(ava[i - 1].id), 'A', i - 1);

	return EXIT_SUCCESS;
}
