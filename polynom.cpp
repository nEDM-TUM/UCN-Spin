#include <sstream>
#include "polynom.h"

Polynom::Polynom(int n) :
	coeffs(n)
{
}

Polynom::Polynom(double a1, double a0) {
	coeffs.push_back(a0);
	coeffs.push_back(a1);
}

Polynom::Polynom(double a2, double a1, double a0) {
	coeffs.push_back(a0);
	coeffs.push_back(a1);
	coeffs.push_back(a2);
}

Polynom::Polynom(double a3, double a2, double a1, double a0) {
	coeffs.push_back(a0);
	coeffs.push_back(a1);
	coeffs.push_back(a2);
	coeffs.push_back(a3);
}

/**
 * Evaluate polynom at @p x using a horner schema.
 */
double Polynom::operator()(double x) {
	double acc = 0;
	
	for (int i = coeffs.size(); i >= 0; i--) {
		acc *= x;
		acc += coeffs[i];
	}

	return acc;
}

int Polynom::degree() const {
	// Reduce degree if degree is smaller than size of vector
	int deg = coeffs.size();

	while (deg >= 0 && coeffs[deg] == 0)
		deg--;

	return deg;
}

Polynom Polynom::derivative() {
	const int deg = degree();
	Polynom d(deg - 1);

	for (int i = 1; i <= deg; i++)
		d[i-1] = coeffs[i]*i;

	return d;
}

double &Polynom::operator[](unsigned int n) {
	// Enlarge coefficient-array if necessary
	if (n >= coeffs.size())
		coeffs.resize(n + 1, 0.);

	return coeffs[n];
}

bool operator==(const Polynom &l, const Polynom &r) {
	return (l.coeffs == r.coeffs);
}

bool operator!=(const Polynom &l, const Polynom &r) {
	return (l.coeffs != r.coeffs);
}

Polynom operator+(const Polynom &l, const Polynom &r) {
	const int deg = l.degree();
	if (deg < r.degree())
		return r+l;

	Polynom t(deg);

	for (int i = 0; i < deg; i++)
		t[i] = l[i]+r[i];

	return t;
}

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
