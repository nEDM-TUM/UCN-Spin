#include "bfield.h"
#include "globals.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include "parameters.h"
#include "basetracking.h"

Bfield::Bfield(const Parameters &theParameters, Basetracking* const btr)
: B0(0.0), B1(0.0), mu0(1.2566371e-6), R(0.0), E0(0.0), g(0.0), ghalf(0.0), B1_g(0.0),
  B1_g_half(0.0), gyroelect(0.0), omegaEDM(0.0), flipangle(0.0), 
  omegalarmor(0.0), factor(-3.9293341e-34*1.5e18), B000(0.0), I(0.0), h2(0.0),
  tracking(btr)
{
	B0 = theParameters.getDoubleParam("B0");
	B1 = theParameters.getDoubleParam("B1");
	E0 = theParameters.getDoubleParam("EfieldMag");
	h2 = theParameters.getDoubleParam("SolenoidHeight") / 2.;
	R = theParameters.getDoubleParam("SolenoidRadius");
	g = theParameters.getDoubleParam("B0Gradient");
	ghalf = g/2.0;
	B1_g = theParameters.getDoubleParam("B1Gradient");
	B1_g_half = B1_g/2.0;
	gyroelect = theParameters.getDoubleParam("GyroelectricRatio");
	omegaEDM = gyroelect*E0;
	flipangle = theParameters.getDoubleParam("Flipangle")/180.0*M_PI;
	omegalarmor = -PRECES_N*B0;
	xyz0[0] = theParameters.getDoubleParam("GradientOffsetX");
	xyz0[2] = theParameters.getDoubleParam("GradientOffsetZ");
	B000 = theParameters.getDoubleParam("SolenoidField");
	I = theParameters.getDoubleParam("SolenoidCurrent");
	#pragma omp master
	{
		cout << "The B0-field = " << B0 << " T" << endl;
		cout << "The B1-field = " << B1 << " T" << endl;
		cout << "The B0-gradient = " << g << " T/m" << endl;
		cout << "The B1-gradient = " << B1_g << " T/m" << endl;
		cout << "The flowing current = " << I << " A" << endl;
		cout << "The Radius of the electrode-circle = " << R << " m" << endl;
		cout << "The heigth of the electrodes = " << h2*2. << " m" << endl;
		cout << "The electric field = " << E0 << " V/m" << endl;
		cout << "The gyroelectric Ratio = " << gyroelect << " m/(V s)" << endl;
		cout << "The flipangle = " << flipangle << " rad = " << flipangle/M_PI*180 << "deg" << endl;
	}
}

Bfield::~Bfield()
{
}

Threevector Bfield::operator()(const double time) const
{
	return eval(time);
}

Threevector Bfield::eval(const double time) const
{
	Threevector position = tracking->getPosition(time);
	debug << "In Bfield: position = " << position.toString() << endl;
	double r = sqrt(position[0]*position[0]+position[1]*position[1]);
	double xr = 0.0, yr = 0.0;
	if(r > 0.0){
		xr = position[0]/r;
		yr = position[1]/r;
	}
	// calculate the B-field of a solenoid at this position
	double BR = Br(r,position[2]);
	double BZ = Bz(r,position[2]);
	return Threevector(BR*xr,BR*yr,BZ);
}

double Bfield::Br(const double r, const double z) const
{
	double c1 = (R+r)*(R+r);
	double c2 = (R-r)*(R-r);
	double a1 = R/sqrt((z+h2)*(z+h2)+c1);
	double a2 = R/sqrt((z-h2)*(z-h2)+c1);
	double k1 = sqrt(((z+h2)*(z+h2)+c2)/((r+h2)*(r+h2)+c1));
	double k2 = sqrt(((z-h2)*(z-h2)+c2)/((r-h2)*(r-h2)+c1));
	return B000*(a1*cel(k1,1,1,-1) - a2*cel(k2,1,1,-1));
}

double Bfield::Bz(const double r, const double z) const
{
	double c1 = (R+r)*(R+r);
	double c2 = (R-r)*(R-r);
	double b1 = (z+h2)/sqrt((z+h2)*(z+h2)+c1);
	double b2 = (z-h2)/sqrt((z-h2)*(z-h2)+c1);
	double g = (R-r)/(R+r);
	double k1 = sqrt(((z+h2)*(z+h2)+c2)/((r+h2)*(r+h2)+c1));
	double k2 = sqrt(((z-h2)*(z-h2)+c2)/((r-h2)*(r-h2)+c1));
	return B000*R/(R+r)*(b1*cel(k1,g*g,1,g) - b2*cel(k2,g*g,1,g));
}

double Bfield::cel(double qqc, double pp, double aa, double bb) const
{
	const double ca = 0.000001; // the desired accuracy is the square of ca
	const double pi02 = 1.5707963268;
	if(qqc == 0.)
	{
		throw "failure in CEL";
	}
	double qc = fabs(qqc), a = aa, b = bb, p = pp, e = qc, em = 1.0, f = 1.0, q = 1.0, g = 1.0;
	if(p > 0.)
	{
		p = sqrt(p);
		b /= p;
	}
	else
	{
		f = qc*qc;
		q = 1.-f;
		g = 1.-p;
		f -= p;
		q *= (b-a*p);
		p = sqrt(f/g);
		a = (a-b)/g;
		b = -q/(g*g*p)+a*p;
	}
	f = a;
	a += b/p;
	g = e/p;
	b = 2*(b+f*g);
	p += g;
	g = em;
	em += qc;
	while(fabs(g-qc) > g*ca)
	{
		qc = sqrt(e);
		qc += qc;
		e = qc*em;
		
		f = a;
		a = a+b/p;
		g = e/p;
		b += f*g;
		b += b;
		p += g;
		g = em;
		em += qc;
	}
	return pi02*(b+a*em)/(em*(em+p));
}
