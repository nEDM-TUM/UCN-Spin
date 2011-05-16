#ifndef BFIELD_H
#define BFIELD_H

#include <fstream>

class Threevector;
class Parameters;
class Basetracking;

class Bfield
{
	public:
		Bfield(const Parameters&, Basetracking* const btr);
		~Bfield();
		double Br(const double r, const double z) const;
		double Bz(const double r, const double z) const;
		double cel(double, double, double, double) const;
		Threevector operator()(const double time) const;
		Threevector eval(const double time) const;

	private:
		double B0, B1, mu0, R, E0, g, ghalf, B1_g, B1_g_half, gyroelect, omegaEDM, flipangle, omegalarmor, factor;
		double xyz0[3];
		double B000,I,h2;
		Basetracking* const tracking;
};
#endif
