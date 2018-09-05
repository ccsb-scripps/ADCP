/*
**  This program calculates average, standard deviation, extreme values,
**  and statistical inefficiency for a given data stream.
**
**  Stats 1.1, Copyright (c) 2005, 2007 Alexei Podtelezhnikov
*/

#define VER "Stats 1.1, Copyright (c) 2005, 2007 Alexei Podtelezhnikov\n"
#define USE "Usage: %s [options] [[-f] filein] [-o fileout]\n\
Options:\n\
 -f infile    default is stdin\n\
 -o outfile   default is stdout\n\
 -i NUMBER    frame size for inefficiency calculations\n\
 -m STRING    input format string, default is \"%%lg\"\n"

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int N = 0;
char *format = "%lg";

void parse_input(void)
{
	int i = 0, j = 0, k, imax = 0, imin = 0;
	double x, y = 0.0, av = 0.0, stdev = 0.0, avn = 0.0, ineff = 0.0;
	double min = 0.0, max = 0.0;

	while ((k = scanf(format, &x)) != EOF) {
		if (k == 0) {
			getchar();	/* eat one character if no match */
			continue;
		}
		if (isnan(x))	/* possibly missing data */
			continue;

		i++;
		av += (x - av) / i;
		stdev += (x * x - stdev) / i;

		if (N) {
			y += x;
			if (i % N == 0) {
				y /= N;
				j++;
				avn += (y - avn) / j;
				ineff += (y * y - ineff) / j;
				y = 0.0;
			}
		}

		if (i == 1 || x > max) {
			max = x;
			imax = i;
		}
		if (i == 1 || x < min) {
			min = x;
			imin = i;
		}
	}

	stdev = stdev - av * av;
	ineff = ineff - avn * avn;
	printf("Count = %d\nMean = %g\nDeviation = %g\n"
	       "Min = %g @ %d\nMax = %g @ %d\n",
	       i, av, sqrt(stdev), min, imin, max, imax);
	if (N)
		printf("Inefficiency = %d * %g / %g = %g\n",
		       N, ineff, stdev, N * ineff / stdev);
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
		case 'i':
			N = atoi(argv[i]);
			break;
		case 'm':
			format = argv[i];
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

	parse_input();

	return EXIT_SUCCESS;
}
