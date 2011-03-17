#include "derivatives.h"
#include "globals.h"
#include <math.h>
#include <iostream>

using namespace std;

/**
 * @class Derivatives
 * The differential equation System for Dopr.
 *
 * @see Dopr
 */

Derivatives::Derivatives(Tracking* TR)
: tracker(TR)
{

}

void Derivatives::operator()(const double t, const double y[], double derivs[])
{
	eval(t, y, derivs);
}

/**
 * Save the value of the derivatives into @p derivs[].
 *
 * @param[in] factor TODO: time?!
 * @param[in] y[] current coordinates of the particle
 * @param[out] derivs[] is set to the derivatives
 */
void Derivatives::eval(const double factor, const double y[], double derivs[])
{
	double B[3] = {0.0,0.0,0.0};
	tracker->getB(factor,B);
	
	/**
	 * The first 3 components of @p derivs are the polarization vector
	 * P and follow the Boch equation.
	 */
	derivs[0] = PRECES_N * (y[1]*B[2] - y[2]*B[1]);
	derivs[1] = PRECES_N * (y[2]*B[0] - y[0]*B[2]);
	derivs[2] = PRECES_N * (y[0]*B[1] - y[1]*B[0]);
}

void Derivatives::evalInitial(const double t, const double y[], double derivs[])
{
	double B[3] = {0.0,0.0,0.0};
	tracker->getInitialB(t,B);
	
	derivs[0] = PRECES_N * (y[1]*B[2] - y[2]*B[1]);
	derivs[1] = PRECES_N * (y[2]*B[0] - y[0]*B[2]);
	derivs[2] = PRECES_N * (y[0]*B[1] - y[1]*B[0]);
}

