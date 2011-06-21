#include "debug.h"
#include "exceptions.h"
#include <iostream>
#include <cmath>

/**
 * @namespace Roots
 * Algorithms for root finding.
 */
namespace Roots {
	template <class T, class C>
	double bisectStep(C* const instance,const T &f, double &x1, double &x2);
	template <class T, class C>
	double safeNewton(C* const instance,const T &f, const T &d, double x1, double x2, double eps);

	const unsigned int MAX_ITERATIONS = 10000;

	class NoRoot : public Exception {
		public:
			NoRoot() : Exception("no root found") {};
	};
}

/**
 * Do one bisection step for <tt>f(double)</tt>. The root must certainly
 * be between @p x1 and @p x2, that means their signs need to be different
 * from each other. At the end of the call, @p x1 or @p x2 will be set to
 * <tt>0.5*(x1+x2)</tt>.
 *
 * @param         f  Function to find roots of
 * @param[in,out] x1 Lower bound of the interval where the root is.
 * @param[in,out] x2 Upper bound of the interval where the root is.
 *
 * @returns The center of the new interval <tt>0.5*(x1+x2)</tt>
 */
template <class T, class C>
double Roots::bisectStep(C* const instance,const T &f, double &x1, double &x2) {
	const double y1 = (instance->*f)(x1);
	const double y2 = (instance->*f)(x2);
	const double xm = .5*(x1+x2);
	const double ym = (instance->*f)(xm);

	debug << "bisecting: (" << x1 << "," << y1 << ") -- (" << xm << "," << ym << ") -- (" << x2 << "," << y2 << ")" << std::endl;

	if (x1 > x2)
		throw "x1 > x2";

	if (y1*ym < 0)
		x2 = xm; // root between xm and x2
	else if (y2*ym < 0)
		x1 = xm; // root between x1 and xm
	else if (y1 == 0)
		x2 = x1;
	else if (y2 == 0)
		x1 = x2;
	else if (ym == 0) {
		x1 = xm;
		x2 = xm;
	}
	else
		throw NoRoot();

	return .5*(x1+x2);
}

/**
 * Find the root of function @p f with derivative @p d. The root must already be
 * bracketed between @p x1 and @p x2 for this to work. This algorithm is quite
 * robust and will use bisection if the newton algorithm does strange things. It
 * not evaluate @p f or @p d outside of <tt>[x1,x2]</tt> except for debugging prints.
 *
 * @param f		  The function to find a root of
 * @param d       The derivative of @p f
 * @param[in] x1  Lower bound of the inverval that brackets the root.
 * @param[in] x2  Upper bound of the inverval that brackets the root.
 * @param[in] eps Accuracy goal, the algorithm will ensure <tt>|f(x_root)| < eps</tt>
 */
template <class T,class C>
double Roots::safeNewton(C* const instance,const T &f, const T &d, double x1, double x2, double eps) {
	if (fabs((instance->*f)(x1)) < eps)
		return x1;
	debug << "Doing initial bisection step:" << std::endl;
	double x = bisectStep(instance,f, x1, x2); // Do a bisection step first to check arguments for sanity
	double y = (instance->*f)(x);
	double last_y = INFINITY;

	unsigned int iteration = 0;

	while (fabs(y) > eps) {
		// Protect against infinite loops
		if(iteration++ > MAX_ITERATIONS)
			throw EndlessLoop();

		const double xn = x - (instance->*f)(x)/(instance->*d)(x);

		debug << "f(" << x << ") = " << y << ", last was: " << last_y << std::endl;

		// Test if result is ok
		if (xn < x1 || xn > x2) {
			// outside of original rage, bisect
			debug << "Bisecting: " << xn << " out of range." << std::endl;
			x = bisectStep(instance,f, x1, x2);
		}
		else if (fabs((instance->*f)(xn)) >= last_y) {
			// function value got bigger, maybe oscillation or divergence
			debug << "Bisecting: |f(" << xn << ")| >= " << last_y << std::endl;
			x = bisectStep(instance,f, x1, x2);
		}
		else {
			// newton seems to go into the right direction, go on
			debug << "Accepting: " << x << " -> " << xn << std::endl;
			x = xn;
		}
		
		last_y = fabs(y);
		y = (instance->*f)(x);
	}

	debug << "Reached accuracy of " << eps << ", f(" << x << ") = " << y << std::endl;

	return x;
}
