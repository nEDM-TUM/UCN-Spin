#ifndef _INTERPOLATIONPOLYNOMIAL_H
#define _INTERPOLATIONPOLYNOMIAL_H

#include <cassert>
#include "polynom.h"
#include "debug.h"

class InterpolationPolynomial : public Polynom
{
	private:
		/**
		 * Power operator for Mathematica generated expressions.
		 */
		static inline double Power(const double base, const int exponent) {
			double ret = 1;
			for (int i = exponent; i >= 1; i--)
				ret *= base;
			return ret;
		};

	public:
		InterpolationPolynomial(double t0, double y0, double d0, double t1, double y1, double d1) : 
			Polynom(3)
		{
			coeffs[0] = (d1*Power(t0,2)*(-(t0*t1) + Power(t1,2)) - Power(t1,3)*y0 + \
			            t0*(Power(t1,2)*(d0*t1 + 3*y0) + t0*(t1*(-(d0*t1) - 3*y1) + \
			            t0*y1)))/(-Power(t1,3) + t0*(t0*(t0 - 3*t1) + 3*Power(t1,2)));
			coeffs[1] = (-(d0*Power(t1,3)) + d1*t0*(-2*Power(t1,2) + t0*(t0 + t1)) + \
			            t0*(2*d0*t0*t1 - d0*Power(t1,2) + t1*(-6*y0 + 6*y1)))/(-Power(t1,3) + \
			            t0*(t0*(t0 - 3*t1) + 3*Power(t1,2)));
			coeffs[2] = (d1*Power(t1,2) + d0*(t0*(-t0 - t1) + 2*Power(t1,2)) + t1*(3*y0 - \
			            3*y1) + t0*(-2*d1*t0 + d1*t1 + 3*y0 - 3*y1))/(-Power(t1,3) + \
			            t0*(t0*(t0 - 3*t1) + 3*Power(t1,2)));
			coeffs[3] = (d1*t0 + d0*(t0 - t1) - d1*t1 - 2*y0 + 2*y1)/(-Power(t1,3) + \
			            t0*(t0*(t0 - 3*t1) + 3*Power(t1,2)));
		};

};

#endif // _INTERPOLATIONPOLYNOMIAL_H
