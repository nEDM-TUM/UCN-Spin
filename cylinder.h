#ifndef CYLINDER_H
#define CYLINDER_H

#include "basegeometry.h"
#include "threevector.h"

class Cylinder : public Basegeometry
{
	public:
		Cylinder(Random *ran, double radius, double height);
		
		bool contains(const Threevector &x) const;
		bool insideHeight(const Threevector &x) const;
		bool insideRadius(const Threevector &x) const;
		
		bool boundsCheck(const Threevector &x);
		void reflect(Threevector &v, const Threevector &x);

		void initialize(Threevector &v, Threevector &x);

		double findIntersection(const double t0, const double t1,
				const Polynom &px, const Polynom &py, const Polynom &pz, double eps);

	private:
		void reflectHeight(Threevector &v);
		void reflectRadius(Threevector &v, const Threevector &x);

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
		bool fReflectHeight;
};

#endif // CYLINDER_H
