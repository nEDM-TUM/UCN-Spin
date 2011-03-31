#include "equationtracker.h"

/**
 * Construct a new @p EquationTracker Object. Do not use it
 * without calling initialize(void) first!
 * @param ran    random number generator
 * @param geo    geometry
 */
EquationTracker::EquationTracker(Random *ran, Basegeometry *geo) : 
	Basetracking(ran, geo), fStepSize(0.01) // TODO: read stepsize as param
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
void EquationTracker::rkStep(const double &h, const double &t, Threevector &x, Threevector &v) {
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
void EquationTracker::makeTrack(const double hmax) {
	// time until which the solution is to be calculated
	const double goal = fTime + hmax;

	// Generate track as far as necessary
	while (fTime < goal) {
		// Do one Runge-Kutta step
		rkStep(fStepSize, fTime, fPos, fVel);

		if (!fGeometry.contains(fPos)) {
			// TODO: backtrack to intersection time, reset fPos and fVel and do smaller step
			
			// Reset fPos and fVel to values before rkStep
			fPos = fTrackpositions.back();
			fVel = fTrackvelocities.back();

			// Find the time when the particle left the volume
			double h = fGeometry.findIntersection(fTime, fTime + fStepSize); // TODO: add interpolated track as param

			// Do smaller step and advance time
			rkStep(h, fTime, fPos, fVel);
			fTime += h;

			// Save velocity and position before collision
			fTracktimes.push_back(fTime);
			fTrackpositions.push_back(fPos);
			fTrackvelocities.push_back(fVel);

			// Reflect
			fGeometry.reflect();

			// Save velocity and position after collision
			fTracktimes.push_back(fTime);
			fTrackpositions.push_back(fPos);
			fTrackvelocities.push_back(fVel);
		}
		else {
			// Step was still inside, just advance time
			fTime += fStepSize;
		}

		// Save new trackpoint
		fTracktimes.push_back(fTime);
		fTrackpositions.push_back(fPos);
		fTrackvelocities.push_back(fVel);
	}
}
