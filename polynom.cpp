#include <sstream>
#include <algorithm>
#include <cassert>
#include "debug.h"
#include "polynom.h"

/**
 * @class Polynom
 *
 * Represents a polynomial.
 */

/**
 * Create a Polynom object ready to take a Polynom of degree @p n without
 * needing a resize of the internal coefficient vector. The coefficients
 * will be set to zero.
 */
Polynom::Polynom(int n) :
	coeffs(n+1, 0.)
{
}

/**
 * Create Polynom <tt>f(x) = a1*x + a0</tt>
 */
Polynom::Polynom(double a1, double a0) :
	coeffs(2)
{
	coeffs[0] = a0;
	coeffs[1] = a1;
}

/**
 * Create Polynom <tt>f(x) = a2*x^2 + a1*x + a0</tt>
 */
Polynom::Polynom(double a2, double a1, double a0) :
	coeffs(3)
{
	coeffs[0] = a0;
	coeffs[1] = a1;
	coeffs[2] = a2;
}

/**
 * Create Polynom <tt>f(x) = a3*x^3 + a2*x^2 + a1*x + a0</tt>
 */
Polynom::Polynom(double a3, double a2, double a1, double a0) :
	coeffs(4)
{
	coeffs[0] = a0;
	coeffs[1] = a1;
	coeffs[2] = a2;
	coeffs[3] = a3;
}

/**
 * Evaluate polynom at @p x using a horner schema.
 */
double Polynom::operator()(double x) const {
	debug << "Polynom: " << this->toString() << " called for x = " << x << std::endl;
	double acc = 0;
	
	for (int i = degree(); i >= 0; i--) {
		debug << "i = " << i << ", acc = " << acc << std::endl;
		acc *= x;
		debug << "acc *= x -> acc = " << acc << std::endl;
		acc += coeffs[i];
		debug << "acc += " << coeffs[i] << " -> acc = " << acc << std::endl;
	}

	debug << "Polynom: result = " << acc << std::endl;
	return acc;
}

/**
 * Return the degree of the polynomial.
 * @returns the biggest n with <tt>coeffs[n] != 0</tt>.
 */
unsigned int Polynom::degree() const {
	if (coeffs.size() == 0)
		return 0;
	assert(coeffs.size() > 0);

	// Reduce degree if degree is smaller than size of vector
	unsigned int deg = coeffs.size() - 1;

	while (deg > 0 && coeffs[deg] == 0)
		deg--;

	assert(deg >= 0);

	return deg;
}

/**
 * Return the derivative.
 */
Polynom Polynom::derivative() const {
	const int deg = degree();
	Polynom d(deg - 1);

	for (int i = 1; i <= deg; i++)
		d[i-1] = coeffs[i]*i;

	return d;
}

/**
 * Get or set coefficients of Polynom.
 * If @p n is greater than the current degree of the
 * polynomial, it is resized.
 */
double &Polynom::operator[](unsigned int n) {
	// Enlarge coefficient-array if necessary
	if (n >= coeffs.size())
		coeffs.resize(n + 1, 0.);

	return coeffs[n];
}

/**
 * Compare two Polynom objects for equality.
 */
bool operator==(const Polynom &l, const Polynom &r) {
	const unsigned int deg = l.degree();

	if (deg != r.degree())
		return false;

	for (int i = deg; i >= 0; i--)
		if (l[i] != r[i])
			return false;

	return true;
}

/**
 * Compare two Polynom objects for inequality.
 */
bool operator!=(const Polynom &l, const Polynom &r) {
	return !(l == r);
}

/**
 * Add two Polynom objects.
 */
Polynom operator+(const Polynom &l, const Polynom &r) {
	// l should be the polynomials with the lower degree, r the one with the higher (or equal) degree. See below for 
	// how this is made sure.

	// common_deg must be the degree that both polynomials have in common, that is the lower of the two degrees
	const unsigned int common_deg = l.degree();
	// max_deg must be the greater of the two degrees
	const unsigned int max_deg = r.degree();

	// Make sure l is the polynomial of lower degree (see above)
	if (common_deg > max_deg)
		return r+l; // Change roles of r and l

	debug << "Adding polynomials " << l.toString() << "  +  " << r.toString() << ", degrees: " << common_deg << ", " << max_deg << std::endl;

	assert(common_deg <= max_deg);

	// The resulting polynomial will have the maximum degree of both
	Polynom t(max_deg);

	for (unsigned int i = 0; i <= max_deg; i++) {\
		if (i <= common_deg) {
			t[i] = l[i]+r[i]; // Add both polynomials
		}
		else {
			assert(i > l.degree());
			assert(i <= r.degree());
			t[i] = r[i]; // Only r[i] has entries of this degree
		}
	}

	return t;
}

/**
 * Substract two polynomials.
 */
Polynom operator-(const Polynom &l, const Polynom &r) {
	return l + (-r);
}

/**
 * Negate a polynomial.
 */
Polynom Polynom::operator-() const {
	Polynom t(degree());
	for (unsigned int i = 0; i < degree(); i++)
		t[i] = -coeffs[i];
	return t;
}

/**
 * Multiply two polynomials
 */
Polynom operator*(const Polynom &l, const Polynom &r) {
	// Get degrees of polynomials
	const int dl = l.degree();
	const int dr = r.degree();

	// Create temporary new polynomial
	Polynom res(dl+dr);
	
	debug << "Created temporary Polynomial: " << res.toString() << std::endl;
	for (unsigned int i = 0; i < res.coeffs.size(); i++) {
		debug << "coeffs[" << i << "] = " << res.coeffs[i] << std::endl;
	}
	debug << "Degree of result polynomial is: " << res.degree() << std::endl;

	// Multiply polynomials
	for (int il = 0; il <= dl; il++)
		for (int ir = 0; ir <= dr; ir++)
			res[il+ir] += r[ir] * l[il];

	debug << "Multiplication returning: " << res.toString() << std::endl;
	return res;
}

/**
 * Return a human-readable string representation.
 */
std::string Polynom::toString() const {
	bool first = true;
	std::ostringstream o;

	for (int i = degree(); i >= 0; i--) {
		const double c = coeffs[i];

		// coefficient
		if (c == 0 && i != 0) {
			continue;
		}
		else if (first) {
			if (c == 1) {
				// do nothing
			}
			else if (c == -1) {
				o << "-";
			}
			else {
				o << c;
			}
			first = false;
		}
		else if (c == 1 && i != 0) {
			o << " + ";
		}
		else if (c < 0) {
			o << " - " << -c;
		}
		else {
			o << " + " << c;
		}
			

		// x^n
		if (i > 1)
			o << "x^" << i;
		else if (i == 1)
			o << "x";
		// do nothing for i == 0
	}

	return o.str();
}
