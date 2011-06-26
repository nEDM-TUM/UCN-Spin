#include <cassert>
#include <limits>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_poly.h>

#include "fastcylindertracker.h"
#include "exceptions.h"

FastCylinderTracker::FastCylinderTracker(const Parameters &params, Random* ran, Cylinder *geo) :
	Basetracking(ran, geo), fGeometry(geo), H(params.getDoubleParam("CylinderHeight")),
	R(params.getDoubleParam("CylinderRadius")), R2(R*R), g(-params.getDoubleParam("GravitationConstant")),
	fLastCollisionSurface(None), fDiffuseProbability(params.getDoubleParam("DiffusionProbability")) {}

FastCylinderTracker::~FastCylinderTracker()
{
}

void FastCylinderTracker::initialize()
{
	fLastCollisionSurface = None;
	Basetracking::initialize();
}

Threevector FastCylinderTracker::getPosition(double t)
{
	static Threevector res; // static for caching
	static double t_cached = std::numeric_limits<double>().quiet_NaN(); // time for which result is cached

	// calculate new result, if the cached one 
	if (t != t_cached) {
		const unsigned int i = findIndex(t);
		const Threevector &x = fTrackpositions[i];
		const Threevector &v = fTrackvelocities[i];
		const double dt = t- fTracktimes[i];
		assert(dt >= 0);

		res[0] = x[0] + v[0]*dt;
		res[1] = x[1] + v[1]*dt;
		res[2] = x[2] + dt*(-0.5*g*dt + v[2]);
		t_cached = t;
	}

	return res;
}

Threevector FastCylinderTracker::getVelocity(double t)
{
	const unsigned int i = findIndex(t);
	const Threevector &v = fTrackvelocities[i];
	const double t0 = fTracktimes[i];
	assert(t >= t0);

	return Threevector(v[0], v[1], v[2] + g*(t-t0));
}

void FastCylinderTracker::makeTrack(double t_start, double h)
{
	// get position and time of last intersection
	Threevector pos = fTrackpositions.back();
	Threevector vel = fTrackvelocities.back();
	double t = fTracktimes.back();

	while (t < t_start + h) {
		debug << "t = " << t << ", x = " << pos << ", v = " << vel << std::endl;
		/**
		 * The coefficients used for the polynomials
		 * are described in Roman's thesis.
		 * The polynomials all have the form
		 * \f[ ax^2 + bx + c \]
		 */
		double a, b, c;

		// collision times
		std::vector< std::pair<double,Surface> > collisions;
		collisions.reserve(6); // allocate enough place for six values

		// DON'T change the order of these!
		// collision at bottom: z(t) = 0
		a = -0.5*g;
		b = g*t + vel[2];
		c = t*(a*t - vel[2]) + pos[2];
		addSolutions(collisions, a, b, c, t, Bottom);
		// collision at top: z(t) - H = 0
		c -= H;
		addSolutions(collisions, a, b, c, t, Top);

		// collision at radius: x(t)*x(t)*y(t)*y(t) - R^2 = 0
		a = vel[0]*vel[0] + vel[1]*vel[1];
		const double temp = 2.*(vel[0]*pos[0] + vel[1]*pos[1]);
		b = -2.*t*a + temp; 
		c = -R*R + t*(t*a - temp) + pos[0]*pos[0] + pos[1]*pos[1];
		addSolutions(collisions, a, b, c, t, Radius);

		// particle has left the volume
		if (collisions.size() <= 0)
			throw LeftVolume();

		debug << "Looking for smallest time in future..." << std::endl;

		// Now find the time of the collision which is smallest time
		// which is greater than the time of the last collision.
		// At the same time, the surface of the last collision is
		// set accordingly.
		fLastCollisionSurface = None;
		double t_coll = std::numeric_limits<double>().max();
		for (unsigned int i = 0; i < collisions.size(); i++) {
			const double candidate = collisions[i].first;
			if (candidate > t && candidate < t_coll) {
				t_coll = candidate;
				fLastCollisionSurface = collisions[i].second;
			}
		}
		assert(t_coll != std::numeric_limits<double>().max()); // was changed
		assert(fLastCollisionSurface != None);
		assert(t_coll > t);
		debug << "Collision at t = " << t_coll << std::endl;

		// now go to the place of the collision
		t = t_coll;
		pos = getPosition(t);
		vel = getVelocity(t);

		// and reflect or scatter
		double rand;
		#pragma omp critical
		{
			rand = fRandomgenerator->uniform();
			assert(rand <= 1 && rand >= 0); 
		}
		if (rand > fDiffuseProbability) {
			assert(fLastCollisionSurface != None);
			if (fLastCollisionSurface == Radius)
				fGeometry->reflectRadius(vel, pos);
			else // Top or Bottom
				fGeometry->reflectHeight(vel);
		}
		else {
			// diffuse scattering
			fGeometry->diffuse(vel, pos);
		}

		// save new trackpoint
		fTracktimes.push_back(t);
		fTrackpositions.push_back(pos);
		fTrackvelocities.push_back(vel);
	}
}

/**
 * Find the soultions of \f$ ax^2 + bx + c = 0 \f$ and add them to the vector @p v.
 * The function will take care that only the absolutely greater value is added, if
 * the surface is the same as the surface of last reflection.
 *
 * @param [in] current the surface on which the reflection takes place
 */
inline void FastCylinderTracker::addSolutions(std::vector< std::pair<double, Surface> > &v, const double a, const double b, const double c,
		const double t0, const Surface current) const
{
	double x1, x2; // used to store results
	debug << "Looking for solutions at " << current << ", last was " << fLastCollisionSurface << std::endl;

	// use GSL to solve quadratic equation#
	// nresults will hold the number of results that were obtained
	int nresults = gsl_poly_solve_quadratic(a, b, c, &x1, &x2);
	debug << "Got " << nresults << " roots: " << x1 << ", " << x2 << std::endl;

	if (current == fLastCollisionSurface) {
		assert(current != None);
		/// If the current surface is the surface of the last reflection, we
		/// need to drop the smaller solution which will be equal to zero, but
		/// only approximately.
		if (nresults == 2) {
			// find the value with the absolute value nearest to zero
			// and make it x2.
			if (fabs(x1-t0) < fabs(x2-t0))
				std::swap(x1, x2);
			
			// add x1 to results
			v.push_back(std::make_pair(x1, current));

			// check x2
			assert(fabs(x2-t0) < 1e-5); // quite lax check
		}
		else {
			// if the equation has less than two solutions, it should
			// have exactly one which is close to zero.
			assert(nresults == 1);
			assert(fabs(x1) < 1e-5); // quite lax check
		}
	}
	else {
		/// If it is another surface, just add all results.
		if (nresults >= 1)
			v.push_back(std::make_pair(x1, current));
		if (nresults == 2)
			v.push_back(std::make_pair(x2, current));
	}
}

unsigned int FastCylinderTracker::findIndex(const double time)
{
	assert(fTrackvelocities.size() == fTrackpositions.size() && fTrackpositions.size() == fTracktimes.size());
	// cache last result
	static unsigned int last_index = -1;

#ifndef NDEBUG
	for (unsigned int j = 1; j < fTracktimes.size(); j++)
		assert(fTracktimes[j-1] < fTracktimes[j]);
#endif

	if (fTracktimes.back() <= time) {
		last_index = fTracktimes.size() - 1;
		return last_index;
	}
	assert(fTracktimes.size() > 1);

	debug << "last_index = " << last_index << ", #trackpoints = " << fTracktimes.size() << std::endl;
	// short circuit if cached result is still valid
	if (last_index >= 0 && fTracktimes[last_index] <= time && fTracktimes[last_index + 1] > time) {
		assert(last_index < fTracktimes.size());
		return last_index;
	}

	// bounds of interval
	unsigned int u = fTracktimes.size() - 1;
	unsigned int l = 0;

	while (u - l != 1) {
		assert(u > 0 && l >= 0);
		assert(u < fTracktimes.size() && l < fTracktimes.size());
		assert(u - l > 1);
		const unsigned int m = (u + l)/2;
		assert(fTracktimes[l] <= fTracktimes[m] && fTracktimes[u] > fTracktimes[m]);
		
		if (fTracktimes[m] > time) {
			u = m;
		}
		else {
			l = m;
		}
		assert(fTracktimes[l] <= time && fTracktimes[u] > time);
	}
	
	debug << "bracketed t = " << time << " between l = " << l << ", t[l] = " << fTracktimes[l]
		  << " and u = " << u << ", t[u] = " << fTracktimes[u] << std::endl;
	assert(u - l == 1);
	assert(fTracktimes[l] <= time && fTracktimes[u] > time);

	debug << "index for t = " << time << " is " << l << std::endl;

	last_index = l;
	assert(l >= 0);
	assert(l < fTracktimes.size());
	return l;
}
