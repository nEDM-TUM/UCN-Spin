#include <iostream>
#include <cmath>

namespace Roots {
	template <class T>
	double bisectStep(T &f, double &x1, double &x2);
	template <class T>
	double safeNewton(T &f, T &d, double x1, double x2, double eps);
}

template <class T>
double Roots::bisectStep(T &f, double &x1, double &x2) {
	const double y1 = f(x1);
	const double y2 = f(x2);
	const double xm = .5*(x1+x2);
	const double ym = f(xm);

	if (x1 > x2)
		throw "x1 > x2";

	if (y1*ym < 0)
		x2 = xm; // root between xm and x2
	else if (y2*ym < 0)
		x1 = xm; // root between x1 and xm
	else
		throw "no root!";

	return .5*(x1+x2);
}

template <class T>
double Roots::safeNewton(T &f, T &d, double x1, double x2, double eps) {
	double x = .5*(x1+x2);
	double y = f(x);
	double last_y = INFINITY;
	
	while (fabs(y) > eps) {
		const double xn = x - f(x)/d(x);

		std::clog << "f(" << x << ") = " << y << ", last was: " << last_y << std::endl;

		// Test if result is ok
		if (xn < x1 || xn > x2) {
			// outside of original rage, bisect
			std::clog << "Bisecting: " << xn << " out of range." << std::endl;
			x = bisectStep(f, x1, x2);
		}
		else if (fabs(f(xn)) >= last_y) {
			// function value got bigger, maybe osciallation or divergence
			std::clog << "Bisecting: |f(" << xn << ")| >= " << last_y << std::endl;
			x = bisectStep(f, x1, x2);
		}
		else {
			// newton seems to go into the right direction, go on
			std::clog << "Accepting: " << x << " -> " << xn << std::endl;
			x = xn;
		}
		
		last_y = fabs(y);
		y = f(x);
	}

	return x;
}
