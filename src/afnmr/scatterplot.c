/* Copyright (c) 1982 Gary Perlman (see Copyright file) */
/*LINTLIBRARY*/

#include <stdio.h>
#include <math.h>

#define	W 12
#define	D  3
#define	MAX_PLOT  100
#define	EPSILON (.000000001)
#define scale(x,min,max,width) ((int) (width*((x)-(min))/((max)-(min)+EPSILON)))
#define	print(x,f) printf ("%10s = %f\n", "x",  x);
	int	elt;
	int	col, row;
	int 	plot[MAX_PLOT][MAX_PLOT];

void scatterplot (float *x, float *y, int n, int plotchar, int height, 
                  int width, int border)
	{
	double	min_x = *x, min_y = *y;
	double	max_x = *x, max_y = *y;
	if (height > MAX_PLOT) height = MAX_PLOT;
	if (width > MAX_PLOT) width = MAX_PLOT;
	for (elt = 0; elt < n; elt++)
		{
		if (x[elt] < min_x) min_x = x[elt];
		if (y[elt] < min_y) min_y = y[elt];
		if (x[elt] > max_x) max_x = x[elt];
		if (y[elt] > max_y) max_y = y[elt];
		}
	for (row = 0; row < height; row++)
		for (col = 0; col < width; col++)
			plot[row][col] = 0;
	for (elt = 0; elt < n; elt++)
		plot[scale(y[elt],min_y,max_y,height)][scale(x[elt],min_x,max_x,width)]++;
	if (border)
		{
		putchar ('|');
		for (col = 0; col < width; col++) putchar ('-');
		putchar ('|');
		if (border > 1) printf ("%g", max_y);
		putchar ('\n');
		}
	for (row = height-1; row >= 0; row--)
		{
		if (border) putchar ('|');
		for (col = 0; col < width; col++)
			if (plot[row][col])
				{
				if (plotchar) putchar (plotchar);
				else if (plot[row][col] >= 10) putchar ('*');
				else putchar (plot[row][col]+'0');
				}
			else putchar (' ');
		if (border) putchar ('|');
		putchar ('\n');
		}
	if (border)
		{
		char	dataline[MAX_PLOT];
		putchar ('|');
		for (col = 0; col < width; col++) putchar ('-');
		putchar ('|');
		if (border > 1) printf ("%g", min_y);
		putchar ('\n');
		if (border > 1)
			{
			sprintf (dataline, "%*.3f", width+2, max_x);
			sprintf (dataline, "%-.3f", min_x);
			for (col = 0; col < width+2; col++)
				if (dataline[col] == NULL) dataline[col] = ' ';
			puts (dataline);
			}
		}
	}
