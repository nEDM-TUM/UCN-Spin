#include "cylinder.h"
#include <math.h>
#include <iostream>

/**
 * @class Cylinder
 * Provides methods to find out if a given vector lies inside a cylinder and
 * to reflect it at the cylinders walls.
 */

/**
 * Create a cylinder.
 * @param radius radius of the cylinder
 * @param height height of the cylinder
 */
Cylinder::Cylinder(double radius, double height)
	: fRadius(radius), fRSquared(radius*radius), fHeight(height)
{
	std::clog << "New cylinder, r = " << fRadius << ", h = " << fHeight << ", r^2 = " << fRSquared << std::endl;
}


/**
 * Check if the cylinder contains the point @p x.
 */
bool Cylinder::contains(double x[])
{
	return insideRadius(x) && insideHeight(x);
}

/**
 * Check if @p x is inside the height of the cylinder.
 *
 * This function does not check if it is inside the radius.
 * @see Cylinder::contains(double[])
 */
bool Cylinder::insideHeight(double x[])
{
	return (x[2] > 0) && (x[2] < fHeight);
}


/**
 * Check if @p x is inside the radius of the cylinder.
 *
 * This function does not check if it is inside the height.
 * @see Cylinder::contains(double[])
 */
bool Cylinder::insideRadius(double x[])
{
	return (x[0]*x[0] + x[1]*x[1]) < fRSquared;
}

/**
 * Reflect v[] if it is outside the cylinder and has been inside
 * since the last reflection.
 *
 * @param[in,out] v veclocity vector, will be changed if necessary.
 * @param[in]     x position of particle, read-only
 * @param[in,out] state will be used to keep track if the particle has been
 *                inside the cylinder since the last reflection. You need to
 *                give the same state for each call of reflect for the same
 *                particle. Initialize to <tt>{false, false}</tt> for new particle!
 *
 * @return true if @p v[] has been changed.
 */
bool Cylinder::reflect(double v[], double x[], bool state[2])
{
	std::clog << "# inside reflect, r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << ", z = " << x[2] << std::endl;
	bool changes = false; // will be set to true if the velocity is changed

	/// Check first, if the particle is inside the height of the cylinder.
	if (!insideHeight(x)) {
		std::clog << "#    z outside of cylinder, state is: " << state[0] << std::endl;
		if (!state[0]) {
			/// Only reflect if particle has been inside the height of the
			/// cylinder since last reflection (so check if <tt>(state[0] ==
			/// false)</tt>).

			// TODO: Add diffusion
			reflectHeight(v);
			changes = true;
		}

		/// Set <tt>state[0] = true</tt> to mark that the particle has been reflected in
		/// z-direction and is not to be reflected again until it has been inside the
		/// cylinder again.
		state[0] = true;
	}
	else {
		/// Set <tt>state[0] = false</tt> to mark hat particle has been inside the height
		/// of the cylinder since the last reflection and must be reflected again of it
		/// leaves it.
		state[0] = false;
		std::clog << "# Reset state[0] of radius" << std::endl;
	}

	/// Now check if the particle is inside the radius of the cylinder.
	if (!insideRadius(x)) {
		std::clog << "#    r outside of cylinder, state is: " << state[1] << std::endl;
		if (!state[1]) {
			/// Only reflect if particle has been inside the radius of the
			/// cylinder since last reflection (so check if <tt>(state[1] ==
			/// false)</tt>).

			// TODO: Add diffusion
			reflectRadius(v, x);
			changes = true;
		}

		/// Set <tt>state[1] = true</tt> to mark that the particle has been reflected in
		/// z-direction and is not to be reflected again until it has been inside the
		/// cylinder again.
		state[1] = true;
	}
	else {
		/// Set <tt>state[1] = false</tt> to mark hat particle has been inside the radius
		/// of the cylinder since the last reflection and must be reflected again of it
		/// leaves it.
		state[1] = false;
		std::clog << "# Reset state[1] of z" << std::endl;
	}

	return changes;
}

void Cylinder::reflectHeight(double v[]) {
	v[2] = -v[2];
}

void Cylinder::reflectRadius(double v[], double x[]) {
	std::clog << "!!! Reflecting x-y, r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << std::endl;
	double n[2]; ///< unitiy vector perpendicular to x[1,2]
	double vn[2]; ///< direction vector of v in x-y-plane
	double v_abs = 0; ///< length of v in x-y-plane
	double sp = 0; ///< scalar product of vn and n

	// fill v_old, v_abs and vn
	for (int i = 0; i <= 2; i++) {
		v_abs += v[i]*v[i];
	}
	v_abs = sqrt(v_abs);
	if (v_abs == 0)
		throw "Cylinder::reflectRadius called with v[0]^2 + v[1]^2 = 0";
	vn[0] = v[0] / v_abs;
	vn[1] = v[1] / v_abs;

	// fill n
	double x_abs = sqrt(x[0]*x[0] + x[1]*x[1]);
	if (x_abs == 0)
		throw "Cylinder::reflectRadius called with x[0]^2 + x[1]^2 = 0";
	n[0] = x[1] / x_abs;
	n[1] = -x[0] / x_abs;

	// Calculate sp
	for (int i = 0; i <= 1; i++) {
		sp += vn[i]*n[i];
	}

	// set v to reflected value
	for (int i = 0; i <= 1; i++) {
		v[i] = -v[i] + 2*v_abs*n[i]*sp;
	}
}
