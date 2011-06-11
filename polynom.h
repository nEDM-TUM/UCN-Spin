#ifndef _POLYNOM_H
#define _POLYNOM_H

#include <vector>
#include <string>

class Polynom {
	public:
		Polynom();
		explicit Polynom(int n);
		Polynom(double a1, double a0);
		Polynom(double a2, double a1, double a0);
		Polynom(double a3, double a2, double a1, double a0);
		unsigned int degree() const;
		Polynom derivative() const;
		std::string toString() const;
		virtual double operator()(double x) const;
		double &operator[](unsigned int n);
		double operator[](unsigned int n) const { return coeffs[n]; };
		Polynom operator-() const;
		friend bool operator==(const Polynom &l, const Polynom& r);
		friend bool operator!=(const Polynom &l, const Polynom& r);
		friend Polynom operator+(const Polynom &l, const Polynom& r);
		friend Polynom operator-(const Polynom &l, const Polynom& r);
		friend Polynom operator*(const Polynom &l, const Polynom& r);
	protected:
		std::vector<double> coeffs;
};

#endif //_POLYNOM_H
