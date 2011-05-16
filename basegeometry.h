#ifndef BASEGEOMETRY_H
#define BASEGEOMETRY_H
 
#include "random.h"
#include "threevector.h"
#include "polynom.h"

/**
 * @class Basegeometry
 * @attention This class keeps internal state, use one instance for one specific particle only!
 */

class Basegeometry
{
	public:
		Basegeometry(Random *ran) : fRandom(ran) {};
		virtual ~Basegeometry();

		/**
		 * Initialize @p x to be inside the cylinder
		 */
		virtual void initialize(Threevector &v, Threevector &x) = 0;

		/**
		 * This method must return if a point x is insinde of the cylinder.
		 * @p x The point which should be checked.
		 */
		virtual bool contains(const Threevector &x) const = 0;

		/**
		 * This method should check, if @p x is inside of the cylinder and
		 * set a reflection state. The reflection state is internal to the
		 * geometry class and should state, at which surface a particle should
		 * be reflected when reflect is called. For example, for a cylinder
		 * the reflection state would contain if the particle is to be reflected
		 * at the top/bottom surface (z-direction) or the round cylinder wall
		 * (r-direction) or both.
		 *
		 * @param x the point for which the bounds check is to be done.
		 * @returns if the particle is outside of the cylinder and should be reflected by calling reflect.
		 * @see reflect
		 */
		virtual bool boundsCheck(const Threevector &x) = 0;

		/**
		 * This method reflects a particle based on the reflection state
		 * which has to be set beforehand by boundsCheck. This method
		 * must always be called after a prior call to boundsCheck has
		 * returned <tt>true</tt>. The parameter @p x may be changed between
		 * both the calls, the decision how to reflect must then be based on
		 * the internal reflection state and not on @p x. The reflection state
		 * has to be reset after the reflection.
		 *
		 * @param v The velocity vector which will be reflected.
		 * @param x The current position of the particle.
		 *
		 * @see boundsCheck
		 */
		virtual void reflect(Threevector &v, const Threevector &x) = 0;

		virtual double findIntersection(const double t0, const double t1,
				const Polynom &px, const Polynom &py, const Polynom &pz, double eps) = 0;

	protected:
		Random *fRandom;
};

#endif
