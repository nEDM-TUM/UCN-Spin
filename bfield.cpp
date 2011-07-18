#include "bfield.h"
#include "globals.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include "parameters.h"
#include "basetracking.h"

Bfield::Bfield(const Parameters &theParameters, Basetracking* const btr)
: mu0(1.2566371e-6), R(0.0), gyromag(0.0), factor(-3.9293341e-34*1.5e18),  B000(0.0), h2(0.0), 
  tracking(btr)
{
	R = theParameters.getDoubleParam("SolenoidRadius");
	centercoil1 = Threevector(0., 0., -R/2.);
	centercoil2 = Threevector(0., 0., R/2);
	gyromag = theParameters.getDoubleParam("GyromagneticRatio");
	B000 = theParameters.getDoubleParam("SolenoidField");
	h2 = theParameters.getDoubleParam("SolenoidHeight") / 2.;
	dipoleposition = Threevector(theParameters.getDoubleParam("Dipolepos0"), theParameters.getDoubleParam("Dipolepos1"), theParameters.getDoubleParam("Dipolepos2"));
	dipole = Threevector(theParameters.getDoubleParam("Dipole0"), theParameters.getDoubleParam("Dipole1"), theParameters.getDoubleParam("Dipole2"));
	earthmagneticfield = theParameters.getIntParam("Earthmagneticfield");
	coilfield = theParameters.getIntParam("Coilfield");
	dipolefield = theParameters.getIntParam("Dipolefield");
	constmagneticfield = theParameters.getIntParam("Constmagneticfield");
	
	#pragma omp master
	{
		cout << "Coilfield: " << coilfield << endl;
		cout << "The Radius of the electrode-circle = " << R << " m" << endl;
		cout << "The heigth of the electrodes = " << h2*2. << " m" << endl;
		cout << "Solenoidfield = " << B000 << " T" << endl;
		cout << "Earthmagneticfield: " << earthmagneticfield << endl;
		cout << "Dipolefield: " << dipolefield << endl;
		cout << "Dipoleposition = " << dipoleposition.toString() << endl;
		cout << "Dipole = " << dipole.toString() << " Tm³/[mu_0]" << endl;
		cout << "Constmagneticfield: " << constmagneticfield << endl;

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
	Threevector field;
	field = Threevector (0.,0.,0.);
	if (earthmagneticfield == 1)
		field = field + evalearthmagneticfield();
	if (coilfield == 1)
		field = field + evalcoil(position + (-1)*centercoil1) + evalcoil(position + (-1)*centercoil2); 
	if (dipolefield == 1)
		field = field + evaldipole(position, dipoleposition, dipole);
	if (constmagneticfield == 1)
		field = field + Threevector(0., 0., 1e-6);
	return field;
}

Threevector Bfield::evalcoil(const Threevector &relposition) const
{
	double r = sqrt(relposition[0]*relposition[0]+relposition[1]*relposition[1]);
	double xr = 0.0, yr = 0.0;
	if(r > 0.0){
		xr = relposition[0]/r;
		yr = relposition[1]/r;
	}
	
	double BR = Br(r,relposition[2]);
	double BZ = Bz(r,relposition[2]);
	return Threevector(BR*xr, BR*yr, BZ);
}

Threevector Bfield::evaldipole(const Threevector &position, const Threevector &positiondipole, const Threevector &magneticdipole) const
{
	Threevector field, relposition, normrelposition;
	double magrelposition, coefficient;
	relposition = position + (-1) * positiondipole;
	magrelposition = relposition.mag();
	normrelposition = relposition.normalized();
	coefficient = mu0 / 4. /3.1459 /magrelposition/magrelposition/magrelposition;
	field = coefficient * (3*(normrelposition*magneticdipole)*normrelposition + (-1)* magneticdipole);
	return field;
}

Threevector Bfield::evalearthmagneticfield () const
{
	return Threevector(17e-6, -30e-6, 3e-6);
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
