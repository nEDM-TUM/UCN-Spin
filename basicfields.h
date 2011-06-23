#ifndef _BASICFIELDS_H
#define _BASICFIELDS_H

#include <string>

#include "bfield.h"
#include "threevector.h"
#include "globals.h"

class DipoleField : public Bfield
{
	public:
		DipoleField(Basetracking *t, const Threevector m, const Threevector pos) : 
			Bfield(t), m(m), dr(-pos) {};
		Threevector eval(const double time) const {
			/// r is set to the vector connecting the position of the dipole x0 with the
			/// current position of the particle x.
			/// \f[ \vec r = \vec x - \vec x_0 \f]
			Threevector r = fTracker->getPosition(time) + dr;
			const double r2 = r.magsquare();
			const double r5 = r2*r2*r.mag();
			// 1e-7 = mu0 / 4pi
			return (1e-7/r5) * (3*r*(m*r)-m*r2);
		};
	private:
		const Threevector m, dr;

};

class HomogenousMagneticField : public Bfield
{
	public:
		HomogenousMagneticField(Basetracking *t, const Threevector B) : 
			Bfield(t), B(B) {};
		Threevector eval(const double time) const {
			return B;
		};
	private:
		const Threevector B;
};

class RelativisticField : public Bfield
{
	public:
		RelativisticField(Basetracking *t, const Threevector E) : 
			Bfield(t), E(E) {};
		Threevector eval(const double time) const {
			const Threevector v = fTracker->getVelocity(time);
			return (-1/(speed_of_light*speed_of_light)) * v.cross(E);
		};
	private:
		const Threevector E;
};

#endif // _BASICFIELDS_H
