int erfcfun_(double *x,double *erf)
{
	double erfc();

	*erf = erfc(*x);
}
