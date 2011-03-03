#ifndef RANDOM_H
#define RANDOM_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class Random
{
	public:
		Random(unsigned long int);
		~Random();
		double uniform();
		double gaussian(const double);
		gsl_rng* GetGsl_Rng();
	private:
		const gsl_rng_type* T;
		gsl_rng* r;
};
#endif
