#ifndef _INTERPOLATIONPOLYNOMIAL_H
#define _INTERPOLATIONPOLYNOMIAL_H

#include <cassert>
#include "polynom.h"
#include "debug.h"

class InterpolationPolynomial : public Polynom
{
	public:
		InterpolationPolynomial(double t0, double y0, double d0, double t1, double y1, double d1) : 
			Polynom(2)
		{
			const double t1_squared = t1*t1;
			const double t0_squared = t0*t0;
			coeffs[0] = (d0*t0*(t0*t1 - t1_squared) + t1_squared*y0 + t0*(-2*t1*y0 + \
			            t0*y1))/(t0*(t0 - 2*t1) + t1_squared);
			coeffs[1] = (d0*(-t0_squared + t1_squared) + t0*(2*y0 - 2*y1))/(t0*(t0 - 2*t1) + t1_squared);
			coeffs[2] = (d0*(t0 - t1) - y0 + y1)/(t0*(t0 - 2*t1) + t1_squared);
		};

};

#endif // _INTERPOLATIONPOLYNOMIAL_H
