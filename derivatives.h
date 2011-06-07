#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "parameters.h"

class Bfield;

class Derivatives
{
	public:
		Derivatives(const Parameters&, Bfield*);
		void operator()(const double, const double[], double[]);
		void eval(const double, const double[], double[]);
	private:
		Bfield *field;
		double gyromag;
};
#endif
