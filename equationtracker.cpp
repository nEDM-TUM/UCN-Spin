#include "equationtracker.h"

/**
 * Construct a new @p EquationTracker Object. Do not use it
 * without calling initialize(void) first!
 * @param ran    random number generator
 * @param geo    geometry
 */
EquationTracker::EquationTracker(Random *ran, Basegeometry *geo) : 
	Basetracking(ran, geo)
{
}

/**
 * Initialize the EquationTracker. Must be called exactly once before any
 * other method of EquationTracker is called.
 */
void EquationTracker::initialize() {
	Basetracking::initialize();
}

/**
 * Do one step of the Runge-Kutta algorithm by
 * integrating derivs.
 *
 * @param[in]     h size of time step that should be taken
 * @param[in]     t current time
 * @param[in,out] x current location, will be updated to new location
 * @param[in,out] v current velocity, will be updated to new velocity
 */
void EquationTracker::rkStep(const double h, const double t, Threevector &x, Threevector &v) {
	static Threevector k1x, k2x, k3x, k4x;
	static Threevector k1v, k2v, k3v, k4v;

	derivs(t, x, v, k1x, k1v);
	derivs(t + h/2., x + .5*k1x, v + .5*k1v, k2x, k2v);
	derivs(t + h/2., x + .5*k2x, v + .5*k2v, k3x, k3v);
	derivs(t + h, x + k3x, v + k3v, k4x, k4v);

	x += (h/6.)*k1x + (h/3.)*k2x + (h/3.)*k3x + (h/6.)*k4x;
	v += (h/6.)*k1v + (h/3.)*k2v + (h/3.)*k3v + (h/6.)*k4v;
}

/**
 * Ensure that track was already generated.
 */
void EquationTracker::makeTrack(const double) {
}
