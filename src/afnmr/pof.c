/* Copyright (c) 1982 Gary Perlman (see Copyright file) */
/*LINTLIBRARY*/

#include <math.h>
double
pof (F, m, n) double F;
	{
	int	i, j;
	int	a, b;
	double	w, y, z, d, p;
	if (F <= 0.0 || m <= 0 || n <= 0) return (1.0);
	a = m%2 ? 1 : 2;
	b = n%2 ? 1 : 2;
	w = F * m / n;
	z = 1.0 / (1.0 + w);
	if (a == 1)
		if (b == 1)
			{
			p = sqrt (w);
			y = 0.3183098862;
			d = y * z / p;
			p = 2 * y * atan (p);
			}
		else
			{
			p = sqrt (w * z);
			d = 0.5 * p * z / w;
			}
	else if (b == 1)
		{
		p = sqrt (z);
		d = 0.5 * z * p;
		p = 1.0 - p;
		}
	else
		{
		d = z * z;
		p = w * z;
		}
	y = 2.0 * w / z;
	for (j = b + 2; j <= n; j += 2)
		{
		d *= (1.0 + a / (j - 2.0)) * z;
		p = (a == 1 ? p + d * y / (j - 1.0) : (p + w) * z);
		}
	y = w * z;
	z = 2.0 / z;
	b = n - 2;
	for (i = a + 2; i <= m; i += 2)
		{
		j = i + b;
		d *= y * j / (i - 2.0);
		p -= z * d / j;
		}
	return (1.0-p);
	}
