#include "bfield.h"
#include "globals.h"
#include <iostream>
#include <math.h>
#include <cstdlib>

Bfield::Bfield(string Bfilename,const Parameters &theParameters)
: B0(0.0), B1(0.0), mu0(1.2566371e-6), R(0.0), E0(0.0), g(0.0), ghalf(0.0), B1_g(0.0),
  B1_g_half(0.0), gyroelect(0.0), omegaEDM(0.0), flipangle(0.0), 
  omegalarmor(0.0), factor(-3.9293341e-34*1.5e18), B000(0.0), I(0.0), h2(0.0),
  Brotationtimes(NULL), NBrottimes(0)
{
	B0 = theParameters.getDoubleParam("B0");
	B1 = theParameters.getDoubleParam("B1");
	E0 = theParameters.getDoubleParam("EfieldMag");
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
	B000 = I*mu0/M_PI/h2/2.;
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
	ifstream Binputfile;
	Binputfile.open(Bfilename.c_str());
	if(!Binputfile)
	{
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}
	double dummy=0.0;
	Binputfile >> dummy;
	while(!(Binputfile.eof()))	
	{
		NBrottimes++;
		Binputfile >> dummy;
	}
	if(NBrottimes%2)
	{
		cerr << "Error: File " << Bfilename << " must contain a number of values that is divisible by 2!" << endl;
		exit(1);
	}
	Binputfile.close();
	ifstream Binputfile2(Bfilename.c_str());
	Brotationtimes = new double[NBrottimes];
	int i=0;
	for(i=0; i<NBrottimes; i++)
	{
		Binputfile2 >> Brotationtimes[i];
	}
	Binputfile2.close();
}

Bfield::~Bfield()
{
	delete[] Brotationtimes;
}

void Bfield::operator()(double t, double *xyz, double *B)
{
	eval(t,xyz,B);
}

void Bfield::eval(double t, double *xyz, double* B)
{
	double co = cos(omegalarmor*t);
	double si = sin(omegalarmor*t);
	double Brot[3] = {0.0,0.0,0.0};
	FieldsRotatingWithLarmorfreq(t, xyz, Brot, co, si);
	NonrotatingFields(xyz, B, co, si);
	B[0] += Brot[0];
	B[1] += Brot[1];
	B[2] += Brot[2];
}

void Bfield::FieldsRotatingWithLarmorfreq(double &t, double *xyz, double *B, double &co, double &si)
{
	B[0] = 0.0;
	B[1] = 0.0;
	B[2] = 0.0;
	if(t <= Brotationtimes[NBrottimes-1])
	{
		int i=0;
		while(t>Brotationtimes[i+1] && i<NBrottimes-1)
		{
			i += 2;
		}
		if((t<Brotationtimes[i+1]) && (t>Brotationtimes[i]))
		{
			if(i==2)
			{
				if(B1_g)
				{
					double xrotated =  co*xyz[0] + si*xyz[1];
					double yrotated = -si*xyz[0] + co*xyz[1];
					B[0] += B1_g*xrotated - B1;
					B[1] += -B1_g_half*yrotated;
					B[2] += -B1_g_half*xyz[2];
				}
				else
				{
					B[0] -= B1;
				}
			}
			else
			{
				if(B1_g)
				{
					double xrotated =  co*xyz[0] + si*xyz[1];
					double yrotated = -si*xyz[0] + co*xyz[1];
					B[0] += B1_g*xrotated + B1;
					B[1] += -B1_g_half*yrotated;
					B[2] += -B1_g_half*xyz[2];				
				}
				else
				{
					B[0] += B1;
				}
			}
		}
	}
	if(E0 && gyroelect)
	{
		B[0] += fabs(gyroelect/PRECES_N)*E0;
	}
}

void Bfield::NonrotatingFields(double *xyz, double *B, double &co, double &si)
{
	double Btemp[3] = {0.0,0.0,0.0};
	if(g)
	{
		Btemp[0] = -ghalf*(xyz[0]-xyz0[0]);
		Btemp[1] = -ghalf*(xyz[1]-xyz0[1]);
		Btemp[2] = g*(xyz[2]-xyz0[2]);
	}
	if(I)
	{
		double BR=0.0, xr=0.0, yr=0.0;
		double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
		if(r > 0.)
		{
			xr = xyz[0]/r;
			yr = xyz[1]/r;
			BR = Br(r,xyz[2]);
		}
		Btemp[0] += xr*BR;
		Btemp[1] += yr*BR;
		Btemp[2] += Bz(r,xyz[2]);
	}
	B[0] =  co*Btemp[0] + si*Btemp[1];
	B[1] = -si*Btemp[0] + co*Btemp[1];
	B[2] = Btemp[2];
}

double Bfield::Br(double r, double z)
{
	double c1 = (R+r)*(R+r);
	double c2 = (R-r)*(R-r);
	double a1 = R/sqrt((z+h2)*(z+h2)+c1);
	double a2 = R/sqrt((z-h2)*(z-h2)+c1);
	double k1 = sqrt(((z+h2)*(z+h2)+c2)/((r+h2)*(r+h2)+c1));
	double k2 = sqrt(((z-h2)*(z-h2)+c2)/((r-h2)*(r-h2)+c1));
	return B000*(a1*cel(k1,1,1,-1) - a2*cel(k2,1,1,-1));
}

double Bfield::Bz(double r, double z)
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

double Bfield::cel(double qqc, double pp, double aa, double bb)
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
