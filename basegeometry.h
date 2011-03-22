#ifndef BASEGEOMETRY_H
#define BASEGEOMETRY_H

#include "threevector.h"

class Basegeometry
{
	public:
		Basegeometry();

		virtual void initialize(Threevector &v, Threevector &x) = 0;
		virtual void contains(Threevector &x) = 0;
		virtual void reflect(Threevector &v, Threevector &x) = 0;

	private:
};

#endif
