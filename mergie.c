/*
** Merger for multiple shuffled MPI-tempered output files.
**
** Copyright (c) 2007 Alexei Podtelezhnikov
*/

#include<stdlib.h>
#include<stdio.h>

int main(int argc, char *argv[])
{
	int i, j, files, *pos;
	char *format, line[83];

	FILE **f;

	format = argv[1];	/* first argument is format */
	files = argc - 2;	/* the rest are file names */

	f = (FILE **) malloc(files * sizeof(FILE *));
	pos = (int *) malloc(files * sizeof(int));

	if (f == NULL || pos == NULL) {
		fprintf(stderr, "Memory Failure\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < files; i++) {
		f[i] = fopen(argv[i + 2], "r");
		if (f[i] == NULL) {
			fprintf(stderr, "Input Failure: %s\n", argv[i + 2]);
			exit(EXIT_FAILURE);
		}
	}

	for (i = 0; i < files; i++) {
		while (fgets(line, sizeof(line), f[i]) != NULL
		       && sscanf(line, format, pos + i) != 1)
			pos[i] = EOF;
	}

	while (1) {
		for (i = 0; i < files; i++)
			if (pos[i] != EOF)
				break;

		for (j = i; i < files; i++)
			if (pos[i] != EOF && pos[i] < pos[j])
				j = i;

		if (j == files)
			break; /* we're done */

		printf(format, pos[j]);
		putchar('\n');
		while (fgets(line, sizeof(line), f[j]) != NULL
		       && sscanf(line, format, pos + j) != 1) {
			printf("%s", line);
			pos[j] = EOF;
		}
	}

	return EXIT_SUCCESS;
}
