#ifndef SPINTRACKING_H
#define SPINTRACKING_H

#include "dopr.h"
#include "derivatives.h"

class Spintracking
{
	public:
		Spintracking(double*, double, double, double, Derivatives*, bool);
		~Spintracking();
		void reset(double*, double, double);
		void precess();
		double getStepperTime();
		double getHdid();
		double dense_out(const int,const double,const double);
		inline int getStepsnottaken() {return stepper->getStepsnottaken();};
	
	private:
		double* P;
		double* dPdt;
		double eps;
		Dopr *stepper;
		Derivatives *derivatives;
		bool dense;
};
#endif
