#include "cylinder.h"
#include "roots.h"
#include "debug.h"
#include "parameters.h"
#include "globals.h"
#include <math.h>
#include <iostream>
#include <cassert>

/**
 * @class Cylinder
 * Provides methods to find out if a given vector lies inside a cylinder and
 * to reflect it at the cylinders walls.
 * @attention This class keeps internal state, use one instance for one specific particle only!
 */

/**
 * Create a cylinder.
 * @param radius radius of the cylinder
 * @param height height of the cylinder
 */
Cylinder::Cylinder(const Parameters &params, Random *ran)
	: Basegeometry(ran), fRadius(params.getDoubleParam("CylinderRadius")),
	fRSquared(fRadius*fRadius), fHeight(params.getDoubleParam("CylinderHeight")),
	fReflectRadius(false), fReflectTop(false), fReflectBottom(false),
	fVelocitySigma(sqrt(params.getDoubleParam("Temperature")*boltzmann/params.getDoubleParam("ParticleMass"))),
	fCutoffSquare(params.getDoubleParam("VelocityCutoff")*params.getDoubleParam("VelocityCutoff")),
	fMinDiffusionAngle(-sin(params.getDoubleParam("MinDiffusionAngle")/180.*M_PI))
{
	assert(fRSquared == fRadius*fRadius);
	debug << "New cylinder, r = " << fRadius << ", h = " << fHeight << ", r^2 = " << fRSquared << std::endl;
}

void Cylinder::initialize(Threevector &v, Threevector &x) {
	// initialize x[0] and x[1] to be inside of radius
	do {
		for (int i = 0; i <= 1; i++)
			x[i] = fRandom->uniform_double(-fRadius, fRadius);
	}
	while (!insideRadius(x));

	assert(x[0]*x[0] + x[1]*x[1] < fRadius*fRadius);
	
	// initialize x[2] to be inside of height
	x[2] = fRandom->uniform_double(0, fHeight);

	// initialize velocity randomly
	do {
		for (int i = 0; i < 3; i++)
			v[i] = fRandom->gaussian(fVelocitySigma);
	} while (fCutoffSquare > 0 && v.magsquare() > fCutoffSquare);

	debug << "initialize: x = " << x.toString() << std::endl;
	debug << "initialize: v = " << v.toString() << std::endl;
}


/**
 * Check if the cylinder contains the point @p x.
 */
bool Cylinder::contains(const Threevector &x) const
{
	return insideRadius(x) && insideHeight(x);
}

/**
 * Check if @p x is inside the height of the cylinder.
 *
 * This function does not check if it is inside the radius.
 * @see Cylinder::contains(double[])
 */
bool Cylinder::insideHeight(const Threevector &x) const
{
	return (x[2] > 0) && (x[2] < fHeight);
}


/**
 * Check if @p x is inside the radius of the cylinder.
 *
 * This function does not check if it is inside the height.
 * @see Cylinder::contains(double[])
 */
bool Cylinder::insideRadius(const Threevector &x) const
{
	assert(fRSquared == fRadius*fRadius);
	return (x[0]*x[0] + x[1]*x[1]) < fRSquared;
}

/**
 * Check if particle is inside of cylinder and set reflection
 * state accordingly. Always call reflect if this method has
 * returned true.
 *
 * @param x The point for which the bounds check is made
 * @returns if the particle is outside the bounds of the geometry
 *
 * @see Basegeometry::boundsCheck
 * @see reflect
 */
bool Cylinder::boundsCheck(const Threevector &x) {
	debug << "boundsCheck for x = " << x.toString() << std::endl;
	/// Set reflection flag if x is outside of height or radius.
	fReflectBottom = false;
	fReflectTop = false;

	if (x[2] < 0) {
		fReflectBottom = true;
	}
	else if (x[2] > fHeight) {
		fReflectTop = true;
	}
	fReflectRadius = !insideRadius(x);

	debug << "fReflectTop = " << fReflectTop << ", fReflectBottom = " << fReflectBottom << ", fReflectRadius = " << fReflectRadius << std::endl;

	/// If any of the flags is set, return true to signal that reflect must be called.
	return (fReflectRadius || fReflectTop || fReflectBottom);
}

/**
 * Reflect v if the previous boundsCheck has returned true.
 *
 * @param[in,out] v veclocity vector, will be changed if necessary.
 * @param[in]     x position of particle, read-only
 *
 * @see Basegeometry::reflect
 */
void Cylinder::reflect(Threevector &v, const Threevector &x)
{
	if (fReflectBottom || fReflectTop) {
		reflectHeight(v);

		// reset state
		fReflectBottom = false;
		fReflectTop = false;
	}

	if (fReflectRadius) {
		reflectRadius(v, x);

		// reset state
		fReflectRadius = false;
	}
}

/**
 * Diffusion method for the EquationTracker which uses the
 * saved state in fReflectRadius, fReflectBottom and fReflectTop.
 */
void Cylinder::diffuse(Threevector &v, const Threevector &pos) {
	/// Diffusion is only implemented at the side wall, in other
	/// cases, reflection is used instead.
	if (fReflectTop || fReflectBottom) {
		reflect(v, pos);
	}
	else {
		assert(fReflectRadius);
		Threevector n(pos);
		n[2] = 0; // n points radially outside the cylinder now
		diffuseAtSurface(v, n);
		// reset reflection state
		fReflectRadius = false;
	}
}

/**
 * Calculate the diffusion of @p v at the surface perpendicular to @p n.
 *
 * @param[in,out] v velocity of the particle, will be used and set to new value
 * @param[in]     n normal vector that is perpendicular to the wall and points to the outside
 */
void Cylinder::diffuseAtSurface(Threevector &v, Threevector n) {
	Threevector v_new; // New velocity
	n.normalize();
	
	do {
		#pragma omp critical
		{
			for (int i = 0; i < 3; i++)
				v_new[i] = fRandom->gaussian(1);
		} 
	} while (v_new.normalized() * n > fMinDiffusionAngle);
	/// The scalar product of @p n and the new velocity is checked to find out
	/// if the diffusion goes into the right direction.

	debug << "diffuse at n = " << n << ", v_new = " << v_new << std::endl;

	// factor which has to be multiplied to vx, vy, vz to conserve velocity
	const double v_fact = sqrt(v.magsquare() / v_new.magsquare());

	assert(fabs(v.magsquare() - (v_fact*v_new).magsquare()) < v.magsquare()*0.0001);
	v = v_fact * v_new;
}

void Cylinder::reflectHeight(Threevector &v) {
	debug << "Before reflectHeight: v = " << v.toString() << std::endl;
	v[2] = -v[2];
	debug << "After reflectHeight: v = " << v.toString() << std::endl;
}

void Cylinder::reflectRadius(Threevector &v, const Threevector &x) {
	debug << "!!! Reflecting x-y, r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << std::endl;
	double n[2]; ///< unitiy vector perpendicular to x[1,2]
	double sp = 0; ///< scalar product of v and n

	// fill n
	double x_abs = sqrt(x[0]*x[0] + x[1]*x[1]);
	if (x_abs == 0)
		throw "Cylinder::reflectRadius called with x[0]^2 + x[1]^2 = 0";
	n[0] = x[1] / x_abs;
	n[1] = -x[0] / x_abs;

	// Calculate sp
	for (int i = 0; i <= 1; i++) {
		sp += v[i]*n[i];
	}

	// set v to reflected value
	for (int i = 0; i <= 1; i++) {
		v[i] = -v[i] + 2*n[i]*sp;
	}
}

double Cylinder::findIntersection(double t0, double t1, const Polynom &px, const Polynom &py, const Polynom &pz, double eps)
{
	assert(eps > 0);
	
	bool try_again = true;
one_more_try:
	double t_height = INFINITY;
	double t_radius = INFINITY;
	debug << "Looking for intersection in [" << t0 << "," << t1 << "]" << std::endl;

	try {
		if (fReflectRadius) {
			// Polynomial for radius: x^2 + y^2 - r^2 == 0
			debug << "Constructing r^2(t) polynomial from px = " << px.toString() << " and py = " << py.toString() << " with offset " << fRSquared << std::endl;
			Polynom rint(px*px + py*py); // Calculate x^2 + y^2
			rint[0] -= fRSquared; // Apply offset
			debug << "Looking for radius intersection, r^2(t) = " << rint.toString() << std::endl;

			// Find intersection time
			t_radius = Roots::safeNewton(rint, rint.derivative(), t0, t1, eps);
		}

		if (fReflectBottom) {
			debug << "Looking for bottom intersection, z(t) = " << pz.toString() << std::endl;
			t_height = Roots::safeNewton(pz, pz.derivative(), t0, t1, eps);
		}

		if (fReflectTop) {
			debug << "Looking for top intersection, z(t) = " << pz.toString() << std::endl;
			Polynom top(pz);
			top[0] -= fHeight; // Apply offset
			t_height = Roots::safeNewton(top, top.derivative(), t0, t1, eps);
		}
	}
	catch (Roots::NoRoot) {
		if (try_again) {
			// Try again once with extended Interval
			try_again = false;
			const double delta = (t1-t0)/2;
			assert(delta > 0);
			t0 -= delta;
			t1 += delta;
			debug << "NO ROOT: Trying again with extended interval [" << t0 << "," << t1 << "]" << std::endl;
			goto one_more_try;
		}
		else {
			throw;
		}
	}

	if (t_radius < t_height)
		return t_radius;
	else
		return t_height;
}
