#include "polynom.h"

class InterpolationPolynomial : public Polynom
{
	private:
		const double h;
		const double t0;

	public:
		InterpolationPolynomial(double t0, double y0, double d0, double t1, double y1, double d1) : 
			Polynom(3), h(t1-t0), t0(t0)
		{
			coeffs[0] = y0;
			coeffs[1] = d0*h;
			coeffs[2] = -2.*d0*h - d1*h - 3.*y0 + 3.*y1;
			coeffs[3] = d0*h + d1*h + 2.*y0 - 2.*y1;
		};

		double operator()(double t) const { return Polynom::operator()(t-t0/h); };
};
