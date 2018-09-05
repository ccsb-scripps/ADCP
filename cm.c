/*
**  This program calculates contact map for a given PDB structure.
**  Parsing options and default values are as follows.
**
**  Coma 1.1, Copyright (c) 2004-2010 Alexei Podtelezhnikov
*/

#define VER "Coma 1.1, Copyright (c) 2004-2010 Alexei Podtelezhnikov\n"
#define USE "Usage: %s [options] [[-f] filein] [-o fileout]\n\
Options:\n\
 -f infile    default is stdin\n\
 -o outfile   default is stdout\n\
 -d 6         cut-off distance, default is 6 angstroms, try zero or negative\n\
 -c A         identifier of the chain to be parsed\n\
 -a CA        atom CA or CB, drop-back default is CA\n\
 -s 10        output symbols, default is 10, try identical symbols\n"

#include"params.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"aadict.h"

typedef double vector[3];

char chain = '\0';
char atom[] = "CA";
char sym[] = "10";
double cutoff = 6.0;
vector ca[2730];
char aa[2730];

int parse_input(void)
{
	char line[83], cpos[5];
	int num, pos, npos;
	double xx, yy, zz;

	num = -1;
	pos = -999;
	while (fgets(line, sizeof(line), stdin) != NULL) {
		/* if chain identifier is given check for it */
		if (chain && line[21] != chain)
			continue;

		/* check for the ATOM and HETATM record */
		if ((line[0] != 'A' || line[1] != 'T') &&
		    (line[0] != 'H' || line[3] != 'A' || line[4] != 'T'))
			continue;

		/* always consider CA as a drop-back case, e.g. for GLY */
		if (line[13] != 'C' || line[14] != 'A' || line[15] != ' ')
			if (line[13] != atom[0] || line[14] != atom[1]
			    || line[15] != ' ')
				continue;

		/* scan position and coordinates */
		if (sscanf(line + 22, "%4[- 0123456789]%*4c%8lf%8lf%8lf",
			   cpos, &xx, &yy, &zz) != 4)
			continue;

		npos = atoi(cpos);
		if (npos != pos) {
			num++;
			if (npos != pos + 1 && pos != -999) {
				ca[num][0] = num * 1e10;
				ca[num][1] = num * 1e10;
				ca[num][2] = num * 1e10;
				aa[num] = '!';
				num++;
			}
		}
		pos = npos;

		ca[num][0] = xx;
		ca[num][1] = yy;
		ca[num][2] = zz;
		aa[num] = aa321(line + 17);

		/* fprintf(stderr, "%d ", pos); */
		if (num > sizeof(ca) / sizeof(vector)) {
			fprintf(stderr, "This file is too big! (%d)\n", num);
			break;
		}
	}
	/* fputc('\n', stderr); */
	return num + 1;
}

void write_contacts(int num)
{
	int i, j;
	char *s;
	vector d;
	double ll;

	s = (char *) malloc((num + 1) * sizeof(char));

	if (s == NULL)
		exit(EXIT_FAILURE);

	for (i = 0; i < num; i++) {
		for (j = 0; j < num; j++) {
			if (i == j) {
				s[j] = aa[j];
				continue;
			}

			d[0] = ca[i][0] - ca[j][0];
			d[1] = ca[i][1] - ca[j][1];
			d[2] = ca[i][2] - ca[j][2];

			ll = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

			if (ll < cutoff * cutoff)
				s[j] = sym[0];
			else if (ll < 1e10)
				s[j] = sym[1];
			else
				s[j] = '!';

		}
		printf("%s\n", s);
	}
	free(s);
}

void write_distances(int num)
{
	int i, j;
	vector d;
	double ll;

	for (i = 0; i < num; putchar('\n'), i++)
		for (j = 0; j < num; j++) {
			if (i == j) {
				printf("  0.0");
				//printf("  %3s", aa123(aa[j]));
				continue;
			}

			d[0] = ca[i][0] - ca[j][0];
			d[1] = ca[i][1] - ca[j][1];
			d[2] = ca[i][2] - ca[j][2];

			ll = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

			printf("%5.1f", (cutoff == 0.0) ?
			       sqrt(ll) : cutoff * cutoff / ll);
		}
}

void write_sparse(int num)
{
	int i, j;
	vector d;
	double ll;

	for (i = 0; i < num; i++)
		for (j = 0; j < num; j++) {
			if (i == j)
				continue;

			d[0] = ca[i][0] - ca[j][0];
			d[1] = ca[i][1] - ca[j][1];
			d[2] = ca[i][2] - ca[j][2];

			ll = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

			if (ll < cutoff * cutoff)
				printf(" %d %d\n", i, j);
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
		case 'a':
			atom[0] = argv[i][0];
			atom[1] = argv[i][1];
			break;
		case 'c':
			chain = argv[i][0];
			break;
		case 'd':
			cutoff = strtod(argv[i], (char **) NULL);
			break;
		case 'f':
			freopen(argv[i], "r", stdin);
			break;
		case 'o':
			freopen(argv[i], "w", stdout);
			break;
		case 's':
			sym[0] = argv[i][0];
			sym[1] = argv[i][1];
			break;
		default:
			fprintf(stderr, VER USE, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char *argv[])
{
	int num;

	read_options(argc, argv);

	/* this will read the coordinates of alpha carbon atoms */
	num = parse_input();

	/* this will output the contact matrix */
	if (sym[0] == sym[1])
		write_sparse(num);
	else if (cutoff > 0.0)
		write_contacts(num);
	else
		write_distances(num);

	return EXIT_SUCCESS;
}
