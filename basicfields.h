#ifndef _BASICFIELDS_H
#define _BASICFIELDS_H

#include <string>
#include <cmath>

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

class GradientField : public Bfield
{
	public:
		GradientField(Basetracking *t, const double gradient) :
			Bfield(t), g(gradient) {};
		Threevector eval(const double time) const {
			const Threevector pos = fTracker->getPosition(time);
			const double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
			const double radial_part = r*g/2.;

			return Threevector(radial_part, radial_part, g*pos[2]);
		};
	private:
		const double g;
};

/**
 * LaserField is a homogenous magnetic field inside a tube in
 * x direction with radius r. Outside the tube, there is no field.
 *
 * @param[in] t      ponter to tracker instance
 * @param[in] B      field inside the tube
 * @param[in] height height of center of tube above z = 0
 * @param[in] radius radius of the tube
 *
 * @return B if inside of tube, 0 else
 */
class LaserField : public Bfield
{
	public:
		LaserField(Basetracking *t, const Threevector B, const double height, const double radius) :
			Bfield(t), B(B), h(height), r2(radius*radius), ZERO(0,0,0) {};
		Threevector eval(const double time) const {
			const Threevector pos = fTracker->getPosition(time);
			const double y = pos[1];
			const double z = pos[2] - h;

			if (y*y+z*z <= r2)
				return B;
			else
				return ZERO;
		};
	private:
		const Threevector B;
		const double h, r2;
		const Threevector ZERO;
};

#endif // _BASICFIELDS_H
