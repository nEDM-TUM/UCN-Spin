#include "random.h"
#include "debug.h"

Random::Random(unsigned long int seed)
: T(gsl_rng_mt19937), r(gsl_rng_alloc(T))
{
	debug << "seed: " << seed << std::endl;
	gsl_rng_set (r, seed);
}

Random::~Random()
{
	gsl_rng_free(r);
}

int Random::uniform_int(int max)
{
	return gsl_rng_uniform_int(r, max+1);
}

double Random::uniform_double(double min, double max)
{
	return gsl_ran_flat(r, min, max);
}

unsigned long int Random::generate_seed()
{
	return gsl_rng_get(r);
}

double Random::uniform()
{
	return gsl_rng_uniform(r);
}

double Random::gaussian(const double sigma)
{
	return gsl_ran_gaussian(r, sigma);
}

double Random::exponential(const double mu)
{
	return gsl_ran_exponential(r, mu);
}

gsl_rng* Random::GetGsl_Rng()
{
	return r;
}
