#include <sstream>
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
	double acc = 0;
	
	for (int i = coeffs.size(); i >= 0; i--) {
		acc *= x;
		acc += coeffs[i];
	}

	return acc;
}

/**
 * Return the degree of the polynomial.
 * @returns the biggest n with <tt>coeffs[n] != 0</tt>.
 */
int Polynom::degree() const {
	// Reduce degree if degree is smaller than size of vector
	int deg = coeffs.size() - 1;

	while (deg >= 0 && coeffs[deg] == 0)
		deg--;

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
	const int deg = l.degree();

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
	const int deg = l.degree();
	if (deg < r.degree())
		return r+l;

	Polynom t(deg);

	for (int i = 0; i < deg; i++)
		t[i] = l[i]+r[i];

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

	// Multiply polynomials
	for (int il = 0; il <= dl; il++)
		for (int ir = 0; ir <= dr; ir++)
			res[il+ir] += r[ir] * l[il];

	return res;
}

/**
 * Return a human-readable string representation.
 */
std::string Polynom::toString() {
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
