#ifndef BFIELD_H
#define BFIELD_H

#include <string>
#include <fstream>
using namespace std;

class Bfield
{
	public:
		Bfield(string, string, string, const double, const double, const double, const double, const double, 
			   const double, const double, const double*, const double, const double, const double, bool);
		~Bfield();
		void FieldsRotatingWithLarmorfreq(double&, double*, double*, double&, double&);
		void NonrotatingFields(double*, double*, double&, double&);
		double Br(double, double);
		double Bz(double, double);
		double cel(double, double, double, double);
		void operator()(double, double*, double*);
		void eval(double, double*, double*);

	private:
		const double B0, B1, mu0, R, E0, g, ghalf, B1_g, B1_g_half, gyroelect, omegaEDM, flipangle, omegalarmor, factor;
		const double *xyz0;
		double B000, I, h2, r0_mag;
		double *Brotationtimes, *Crosstalkrotationtimes, *r0;
		int NBrottimes, NCrossrottimes;
		bool rotatedipole;
};
#endif
