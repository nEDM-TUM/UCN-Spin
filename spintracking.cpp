#include "spintracking.h"

using namespace std;

Spintracking::Spintracking(double* PolVec, double initialT, double error, double firsthtry, Derivatives *DERI, bool dens)
: P(PolVec), dPdt(NULL), eps(error), stepper(NULL), derivatives(DERI), dense(dens)
{
	dPdt = new double[3];
	stepper = new Dopr(firsthtry, 3, P, dPdt, initialT, derivatives, eps, eps, dense);
}

Spintracking::~Spintracking()
{
	delete stepper;
	delete[] dPdt;
}

double Spintracking::getStepperTime()
{
	return stepper->getT();
}

double Spintracking::getHdid()
{
	return stepper->getHdid();
}
void Spintracking::reset(double* PolVec, double initialT, double firsthtry)
{
	P = PolVec;
	
	derivatives->evalInitial(initialT, P, dPdt);
	stepper->reset(firsthtry, P, dPdt, initialT);
}

void Spintracking::precess()
{
	stepper->step();
}

double Spintracking::dense_out(const int i,const double t,const double h)
{
	return stepper->dense_out(i,t,h);
}
