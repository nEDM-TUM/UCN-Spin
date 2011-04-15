#include "equationtracker.h"
#include "debug.h"
#include "interpolationpolynomial.h"

/**
 * Construct a new @p EquationTracker Object. Do not use it
 * without calling initialize(void) first!
 * @param ran    random number generator
 * @param geo    geometry
 */
EquationTracker::EquationTracker(Random *ran, Basegeometry *geo) : 
	Basetracking(ran, geo), fStepSize(0.01), fTime(0) // TODO: read stepsize as param
{
}

/**
 * Initialize the EquationTracker. Must be called exactly once before any
 * other method of EquationTracker is called.
 */
void EquationTracker::initialize() {
	fTime = 0;
	Basetracking::initialize();
	fPos = fTrackpositions.back();
	fVel = fTrackvelocities.back();
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
	Threevector k1x, k2x, k3x, k4x;
	Threevector k1v, k2v, k3v, k4v;

	debug << "i\tk1x\tk1v" << std::endl;
	derivs(t, x, v, k1x, k1v);
	k1x *= h;
	k1v *= h;
	debug << "1\t" << k1x.toString() << "\t" << k1v.toString() << std::endl;
	derivs(t + h/2., x + .5*k1x, v + .5*k1v, k2x, k2v);
	k2x *= h;
	k2v *= h;
	debug << "2\t" << k2x.toString() << "\t" << k2v.toString() << std::endl;
	derivs(t + h/2., x + .5*k2x, v + .5*k2v, k3x, k3v);
	k3x *= h;
	k3v *= h;
	debug << "3\t" << k3x.toString() << "\t" << k3v.toString() << std::endl;
	derivs(t + h, x + k3x, v + k3v, k4x, k4v);
	k4x *= h;
	k4v *= h;
	debug << "4\t" << k4x.toString() << "\t" << k4v.toString() << std::endl;

	x += (1./6.)*k1x + (1./3.)*k2x + (1./3.)*k3x + (1./6.)*k4x;
	v += (1./6.)*k1v + (1./3.)*k2v + (1./3.)*k3v + (1./6.)*k4v;
}

/**
 * Ensure that track was already generated.
 */
void EquationTracker::makeTrack(const double hmax) {
	// time until which the solution is to be calculated
	const double goal = fTime + hmax;

	// Generate track as far as necessary
	while (fTime < goal) {
		debug << "TIME: " << fTime << ", dt = " << fStepSize << std::endl;
		debug << "Before step: " << fPos.toString() << " speed " << fVel.toString() << std::endl;
		// Do one Runge-Kutta step
		rkStep(fStepSize, fTime, fPos, fVel);

		debug << "After step: " << fPos.toString() << " speed " << fVel.toString() << std::endl;
		// Check if the particle has left the volume during the step
		if (fGeometry->boundsCheck(fPos)) {
			/// If the particle left the volume during the step, backtrack to the intersection time,
			/// do the reflection and save position and velocity before and after the collision.
			debug << "Particle left Volume, x = " << fPos.toString() << std::endl;

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
			// TODO: make eps flexible
			double t_coll = fGeometry->findIntersection(fTime, t1, px, py, pz, 1e-10);

			// Do smaller step and advance time
			rkStep(t_coll - fTime, fTime, fPos, fVel);
			fTime = t_coll;

			debug << "AT COLLISION: " << fPos.toString() << " speed " << fVel.toString() << std::endl;
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
		debug << "At end of step: " << fPos.toString() << " speed " << fVel.toString() << std::endl;
	}
}

Threevector EquationTracker::getPosition(double time) {
	// TODO
	return Threevector();
}
