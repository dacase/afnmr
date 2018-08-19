/* Copyright (c) 1982 Gary Perlman (see Copyright file) */
/*LINTLIBRARY*/

/* strings reads from ioptr into abase, an array of most maxstrings strings,
   at most maxstrings strings, each of length at most maxchars-1 chars.
   It returns the number of strings read in, or maxstrings + 1 if some
   information is discarded.
*/
#include <stdio.h>
#include <ctype.h>
#ifndef ESCAPE
#define	ESCAPE '\\'
#endif
int sstrings (char *line, char *abase, int maxstrings, int maxchars)
	{
	int	nstrings = 0;
	int	nchars;

	while (isspace (*line)) line++;
	while (*line)
		{
		nchars = 0;
		while (*line && !isspace (*line) && nchars<maxchars-1)
			if ((abase[nchars++] = *line++) == ESCAPE)
				abase[nchars-1] = *line++;
		abase[nchars] = '\0';
		abase += maxchars;
		while (*line && !isspace (*line)) line++;
		while (isspace (*line)) line++;
		if (++nstrings == maxstrings)
			return (maxstrings + (*line ? 1 : 0));
		}
	return (nstrings);
	}

int fstrings (FILE *ioptr, char *abase, int maxstrings, int maxchars)
	{
	char	line[BUFSIZ];

	if (fgets (line, BUFSIZ, ioptr) == NULL) return (EOF);
	return (sstrings (line, abase, maxstrings, maxchars));
	}

