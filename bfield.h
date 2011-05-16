#ifndef BFIELD_H
#define BFIELD_H

#include <string>
#include <fstream>
#include "parameters.h"

using namespace std;

class Bfield
{
	public:
		Bfield(string,const Parameters&);
		~Bfield();
		void FieldsRotatingWithLarmorfreq(double&, double*, double*, double&, double&);
		void NonrotatingFields(double*, double*, double&, double&);
		double Br(double, double);
		double Bz(double, double);
		double cel(double, double, double, double);
		void operator()(double, double*, double*);
		void eval(double, double*, double*);

	private:
		double B0, B1, mu0, R, E0, g, ghalf, B1_g, B1_g_half, gyroelect, omegaEDM, flipangle, omegalarmor, factor;
		double xyz0[3];
		double B000,I,h2;
		double *Brotationtimes;
		int NBrottimes;
};
#endif
