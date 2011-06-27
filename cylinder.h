#ifndef CYLINDER_H
#define CYLINDER_H

#include "basegeometry.h"
#include "threevector.h"
#include "parameters.h"

class FastCylinderTracker;

class Cylinder : public Basegeometry
{
	public:
		Cylinder(const Parameters& params, Random *ran);
		
		bool contains(const Threevector &x) const;
		bool insideHeight(const Threevector &x) const;
		bool insideRadius(const Threevector &x) const;
		
		bool boundsCheck(const Threevector &x);
		void reflect(Threevector &v, const Threevector &x);
		void diffuse(Threevector &v, const Threevector &x);
		void diffuseAtSurface(Threevector &v, Threevector n);

		void initialize(Threevector &v, Threevector &x);

		double findIntersection(double t0, double t1,
				const Polynom &px, const Polynom &py, const Polynom &pz, double eps);

		void reflectHeight(Threevector &v);
		void reflectRadius(Threevector &v, const Threevector &x);

		friend class FastCylinderTracker;

	private:
		double fRadius; ///< radius of cylinder
		double fRSquared; ///< squared radius of cylinder
		double fHeight; ///< height of cylinder
		
		/**
		 * Reflection state for reflect and boundsCheck.
		 * @see Basegeometry::boundsCheck
		 * @see Basegeometry::reflect
		 */
		bool fReflectRadius;
		/**
		 * Reflection state for reflect and boundsCheck.
		 * @see Basegeometry::boundsCheck
		 * @see Basegeometry::reflect
		 */
		bool fReflectTop;
		bool fReflectBottom;

		double fVelocitySigma; ///< sigma for maxwell distribution of velocity
		double fCutoffSquare; ///< highest possible velocity
		double fMinDiffusionAngle; //< minimal angle after diffusion versus surface
};

#endif // CYLINDER_H
