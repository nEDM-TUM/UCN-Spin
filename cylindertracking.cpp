/**
 * @file
 *
 * Make a track for a specific velocity and starting point
 */

#include <iostream>
#include <cmath>
#include "basetracking.h"
#include "cylinder.h"
#include "parameters.h"

/**
 * Find the time of the next intersection with the
 * cylinder walls.
 */
double nextIntersection(double x0[], double v0[], double t0) {
	// times of all intersections that where found
	double times[6];

	// set all the times to infinity as initialisation
	for (int i = 0; i < 6; i++)
		times[i] = INFINITY;

	// squared quantities
	double vx0_squared = v0[0]*v0[0];
	double vy0_squared = v0[1]*v0[1];
	double vz0_squared = v0[2]*v0[2];

	// determinats of quadratic equations
	double detR = R_squared*(vx0_squared+vy0_squared) - (v0[1]*x0[0] - v0[0]*x0[1])*(v0[1]*x0[0] - v0[0]*x0[1]);
	double detB = vz0_squared + 2.*g*x0[2]; // bottom of cylinder
	double detH = vz0_squared + 2.*g*(x0[2] - H); // top of cylinder

	// TODO: calculate times
	
	// TODO: cancel out starting point if track started at wall
	
	// Return minimum of valid times
	double t_min = INFINITY;
	for (int i = 0; i < 6; i++)
		if (times[i] < t_min)
			t_min = times[i];

	return t_min;
}
