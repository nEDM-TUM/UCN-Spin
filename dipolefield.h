#ifndef DIPOLEFIELD_H
#define DIPOLEFIELD_H

class Dipolefield
{
	public:
		Dipolefield(double);
		void getField(double, double *, double *, double *);
		
	private:
		const double factor;
		double *r;
		double omegalarmor;
};

#endif
