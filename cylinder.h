#ifndef CYLINDER_H
#define CYLINDER_H

#include "threevector.h"

class Cylinder
{
	public:
		Cylinder(double radius, double height);
		
		bool contains(const Threevector &x);
		bool insideHeight(const Threevector &x);
		bool insideRadius(const Threevector &x);

		bool reflect(Threevector &v, Threevector &x, bool state[2]);

	private:
		void reflectHeight(Threevector &v);
		void reflectRadius(Threevector &v, Threevector &x);

		double fRadius; ///< radius of cylinder
		double fRSquared; ///< squared radius of cylinder
		double fHeight; ///< height of cylinder
};

#endif // CYLINDER_H
