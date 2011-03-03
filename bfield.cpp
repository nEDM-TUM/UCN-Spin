#include "bfield.h"
#include "globals.h"
#include <iostream>
#include <math.h>
#include <cstdlib>

Bfield::Bfield(string Bfilename, string Crossfilename, string Crosslocfilename, const double B0_value, const double B1_value,
			   const double I_value, const double R_value, const double h_value, const double E0_value, const double g_value, 
			   const double *xyz0_value, const double B1_g_value, const double gyroelect_value, const double flipangle_value, 
			   bool rotatedipole_value)
: B0(B0_value), B1(B1_value), mu0(1.2566371e-6), R(R_value), E0(E0_value), g(g_value), ghalf(g_value/2.), B1_g(B1_g_value),
  B1_g_half(B1_g_value/2.0), gyroelect(gyroelect_value), omegaEDM(gyroelect*E0_value), flipangle(flipangle_value), 
  omegalarmor(-PRECES_N*B0_value), factor(-3.9293341e-34*1.5e18), xyz0(xyz0_value), B000(0.0), I(I_value), h2(h_value/2.), 
  r0_mag(0.0), Brotationtimes(NULL), Crosstalkrotationtimes(NULL), r0(NULL), NBrottimes(0), NCrossrottimes(0), 
  rotatedipole(rotatedipole_value)
{
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
	B000 = I*mu0/M_PI/h2/2.;
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
	
	ifstream Crossinputfile;
	Crossinputfile.open(Crossfilename.c_str());
	if(!Crossinputfile)
	{
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}
	dummy=0.0;
	Crossinputfile >> dummy;
	while(!(Crossinputfile.eof()))	
	{
		NCrossrottimes++;
		Crossinputfile >> dummy;
	}
	if(NCrossrottimes%2)
	{
		cerr << "Error: File " << Crossfilename << " must contain a number of values that is divisible by 2!" << endl;
		exit(1);
	}
	Crossinputfile.close();
	ifstream Crossinputfile2(Crossfilename.c_str());
	Crosstalkrotationtimes = new double[NCrossrottimes];
	i=0;
	for(i=0; i<NCrossrottimes; i++)
	{
		Crossinputfile2 >> Crosstalkrotationtimes[i];
	}
	Crossinputfile2.close();
	
	r0 = new double[3];
	ifstream Crosslocations;
	Crosslocations.open(Crosslocfilename.c_str());
	for(i=0; i<3; i++)
	{
		Crosslocations >> r0[i];
		r0_mag += r0[i]*r0[i];
	}
	r0_mag = sqrt(r0_mag);
	cout << "r0_mag: " << r0_mag << endl;
	Crosslocations.close();
}

Bfield::~Bfield()
{
	delete[] Brotationtimes;
	delete[] Crosstalkrotationtimes;
	delete[] r0;
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
	if(r0_mag)
	{
		if(t <= Crosstalkrotationtimes[NCrossrottimes-1])
		{
			int i=0;
			while(t>Crosstalkrotationtimes[i+1] && i<NCrossrottimes-1)
			{
				i += 2;
			}
			if((t<Crosstalkrotationtimes[i+1]) && (t>Crosstalkrotationtimes[i]))
			{
				double angle = flipangle;
				if(rotatedipole)
				{
					angle += omegaEDM*t;
				}
				double P[3] = {0.0,sin(-angle),cos(angle)};
				double r[3] = {xyz[0]-r0[0],xyz[1]-r0[1],xyz[2]-r0[2]};
				double rrotated[2] = {co*r[0]+si*r[1],-si*r[0]+co*r[1]};
				double r_mag = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
				rrotated[0] /= r_mag;
				rrotated[1] /= r_mag;
				r[2] /= r_mag;
				double DistFact = factor/(r_mag*r_mag*r_mag);
				double ScalProd = 3*(P[0]*rrotated[0]+P[1]*rrotated[1]+P[2]*r[2]);
				B[0] += DistFact*(ScalProd*rrotated[0]-P[0]);
				B[1] += DistFact*(ScalProd*rrotated[1]-P[1]);
				B[2] += DistFact*(ScalProd*r[2]-P[2]);
			}
		}
	}
	if(E0 && gyroelect && !(rotatedipole))
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
