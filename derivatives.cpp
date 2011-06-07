#include "derivatives.h"
#include "globals.h"
#include <math.h>
#include <iostream>
#include "bfield.h"
#include "threevector.h"
#include "parameters.h"

using namespace std;

/**
 * @class Derivatives
 * The differential equation System for Dopr.
 *
 * @see Dopr
 */

Derivatives::Derivatives(const Parameters& theParameters, Bfield* F)
: field(F), gyromag(0.0)
{
	gyromag = theParameters.getDoubleParam("GyromagneticRatio");
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
void Derivatives::eval(const double time, const double y[], double derivs[])
{
	Threevector B = field->eval(time);
	debug << "In eval: B = " << B.toString() << endl;
	
	/**
	 * The first 3 components of @p derivs are the polarization vector
	 * P and follow the Boch equation.
	 */
	derivs[0] = gyromag * (y[1]*B[2] - y[2]*B[1]);
	derivs[1] = gyromag * (y[2]*B[0] - y[0]*B[2]);
	derivs[2] = gyromag * (y[0]*B[1] - y[1]*B[0]);
}

