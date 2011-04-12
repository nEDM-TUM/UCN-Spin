#include "equationtracker.h"
#include "interpolationpolynomial.h"

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

		// Check if the particle has left the volume during the step
		if (fGeometry->boundsCheck(fPos)) {
			/// If the particle left the volume during the step, backtrack to the intersection time,
			/// do the reflection and save position and velocity before and after the collision.

			// Save current vectors and time for interpolation polynomial
			Threevector pos1(fPos);
			Threevector vel1(fVel);
			double t1 = fTime + fStepSize;
			
			// Reset fPos and fVel to values before rkStep
			fPos = fTrackpositions.back();
			fVel = fTrackvelocities.back();

			// Create interpolation polynomials for x, y and z
			InterpolationPolynomial px(fTime, fPos[0], fVel[0], t1, pos1[0], vel1[0]);
			InterpolationPolynomial py(fTime, fPos[1], fVel[1], t1, pos1[1], vel1[1]);
			InterpolationPolynomial pz(fTime, fPos[2], fVel[2], t1, pos1[2], vel1[2]);

			// Find the time when the particle left the volume
			double t_coll = fGeometry->findIntersection(fTime, t1, px, py, pz, fStepSize * fVel.mag());

			// Do smaller step and advance time
			rkStep(t_coll - fTime, fTime, fPos, fVel);
			fTime = t_coll;

			// Save velocity and position before collision
			fTracktimes.push_back(fTime);
			fTrackpositions.push_back(fPos);
			fTrackvelocities.push_back(fVel);

			// Reflect
			fGeometry->reflect(fVel, fPos);
		}
		else {
			// Step was still inside, just advance time
			fTime += fStepSize;
		}

		// Save new trackpoint (in case of collision the one after the refelction)
		fTracktimes.push_back(fTime);
		fTrackpositions.push_back(fPos);
		fTrackvelocities.push_back(fVel);
	}
}
