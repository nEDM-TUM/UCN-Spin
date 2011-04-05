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
		int degree() const;
		Polynom derivative();
		std::string toString();
		double operator()(double x);
		double &operator[](unsigned int n);
		double operator[](unsigned int n) const { return coeffs[n]; };
		Polynom operator-();
		friend bool operator==(const Polynom &l, const Polynom& r);
		friend bool operator!=(const Polynom &l, const Polynom& r);
		friend Polynom operator+(const Polynom &l, const Polynom& r);
		friend Polynom operator-(const Polynom &l, const Polynom& r);
	private:
		std::vector<double> coeffs;
};

#endif //_POLYNOM_H
