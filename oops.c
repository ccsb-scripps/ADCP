/*
**  This program creates a PostScript path of protein structure
**  provided in a PDB file.
**
**  Oops 1.0, Copyright (c) 2004, 2007 Alexei Podtelezhnikov
*/

#define VER "Oops 1.0, Copyright (c) 2004, 2005 Alexei Podtelezhnikov\n"
#define USE "Usage: %s [options] [[-f] filein] [-o fileout]\n\
Options:\n\
  -c A     identifier for a chain to be displayed\n\
  -m MASK  hexadecimal mask for atoms ...-0-C-Cb-Ca-H-N\n"

#include<stdio.h>
#include<stdlib.h>

typedef double vector[3];

char chain = '\0';
int req = 0x4;

void output(int flag, vector nl, vector hl, vector cal, vector cbl, vector cl,
	    vector ol)
{
	if (req & 0x1 && flag & 0x1) {
		printf("[ %g %g %g ]\n", nl[0], nl[1], nl[2]);
		if (req & 0x2 && flag & 0x2) {
			printf("[ %g %g %g ]\n", hl[0], hl[1], hl[2]);
			printf("[ %g %g %g ]\n", nl[0], nl[1], nl[2]);
		}
	}

	if (req & 0x4 && flag & 0x4) {
		printf("[ %g %g %g ]\n", cal[0], cal[1], cal[2]);
		if (req & 0x8 && flag & 0x8) {
			printf("[ %g %g %g ]\n", cbl[0], cbl[1], cbl[2]);
			printf("[ %g %g %g ]\n", cal[0], cal[1], cal[2]);
		}
	}

	if (req & 0x10 && flag & 0x10) {
		printf("[ %g %g %g ]\n", cl[0], cl[1], cl[2]);
		if (req & 0x20 && flag & 0x20) {
			printf("[ %g %g %g ]\n", ol[0], ol[1], ol[2]);
			printf("[ %g %g %g ]\n", cl[0], cl[1], cl[2]);
		}
	}
}

void pdbparse(void)
{
	char line[81];
	char pos[] = "    ";
	int mask, flag = 0;
	double *atm;
	vector nl, hl, cal, cbl, cl, ol;

	while (fgets(line, sizeof(line), stdin) != NULL) {
		/* check for an ATOM record */
		if (line[0] != 'A' || line[1] != 'T')
			continue;

		/* if chain identifier is given check for it */
		if (chain && line[21] != chain)
			continue;

		if (line[25] != pos[3] || line[24] != pos[2] ||
		    line[23] != pos[1] || line[22] != pos[0]) {
			pos[0] = line[22];
			pos[1] = line[23];
			pos[2] = line[24];
			pos[3] = line[25];

			output(flag, nl, hl, cal, cbl, cl, ol);
			flag = 0;
		}

		/* check the atom type */
		if (req & 0x1 && line[13] == 'N' && line[14] == ' ') {
			atm = nl;
			mask = 0x1;
		} else if (req & 0x2 && line[13] == 'H' && line[14] == ' ') {
			atm = hl;
			mask = 0x2;
		} else if (req & 0x4 && line[13] == 'C' && line[14] == 'A') {
			atm = cal;
			mask = 0x4;
		} else if (req & 0x8 && line[13] == 'C' && line[14] == 'B') {
			atm = cbl;
			mask = 0x8;
		} else if (req & 0x10 && line[13] == 'C' && line[14] == ' ') {
			atm = cl;
			mask = 0x10;
		} else if (req & 0x20 && line[13] == 'O' && line[14] == ' ') {
			atm = ol;
			mask = 0x20;
		} else
			continue;

		/* scan coordinates and check for a good scan */
		if (sscanf(line + 30, "%8lf%8lf%8lf", atm, atm + 1, atm + 2)
		    != 3)
			continue;

		flag |= mask;
	}

	output(flag, nl, hl, cal, cbl, cl, ol);
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
		case 'c':
			chain = argv[i][0];
			break;
		case 'm':
			req = strtol(argv[i], (char **) NULL, 16);
			break;
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
	read_options(argc, argv);

	pdbparse();

	return EXIT_SUCCESS;
}
