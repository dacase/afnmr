int erfcfun_(x,erf)
	double *x,*erf;
{
	double erfc();

	*erf = erfc(*x);
}
