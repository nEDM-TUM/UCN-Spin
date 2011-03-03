#include "random.h"

Random::Random(unsigned long int seed)
: T(gsl_rng_mt19937), r(gsl_rng_alloc(T))
{
	gsl_rng_set (r, seed);
}

Random::~Random()
{
	gsl_rng_free(r);
}

double Random::uniform()
{
	return gsl_rng_uniform(r);
}

double Random::gaussian(const double sigma)
{
	return gsl_ran_gaussian(r, sigma);
}

gsl_rng* Random::GetGsl_Rng()
{
	return r;
}
