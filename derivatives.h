#ifndef DERIVATIVES_H
#define DERIVATIVES_H

class Bfield;

class Derivatives
{
	public:
		Derivatives(Basetracking*);
		void operator()(const double, const double[], double[]);
		void eval(const double, const double[], double[]);
	private:
		Bfield *field;
};
#endif
