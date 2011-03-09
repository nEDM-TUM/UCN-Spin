#include "tracking.h"
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;


Tracking::Tracking(Random *ran, Bfield* BF, Parameters& theParameters)
: B(BF), rand(ran), Pos(NULL), xyz1(NULL), xyz2(NULL), dir1(NULL), dir2(NULL), posxyz(NULL), sqrt2D(0.00025), 
  lambda(0.0), Rsquare(0.0), R(0.0), scaling(1.), t(0.0), usetwopoints(false), htry(0.0)
{
	R = theParameters.getDoubleParam("Radius");
	Rsquare = R*R;
	Pos = new double[3];
	xyz1 = new double[3];
	xyz2 = new double[3];
	dir1 = new double[3];
	dir2 = new double[3];
	posxyz = new double[3];
	for(int j=0; j<3; j++)
	{
		Pos[j] = 0.0;
		xyz1[j] = 0.0;
		xyz2[j] = 0.0;
		dir1[j] = 0.0;
		dir2[j] = 0.0;
		posxyz[j] = 0.0;
	}
}

Tracking::~Tracking()
{
	delete[] Pos;
	delete[] xyz1;
	delete[] xyz2;
	delete[] dir1;
	delete[] dir2;
	delete[] posxyz;
}

void Tracking::getB(double factor, double *Bvec)
{
	double factofwholestep = factor*scaling;
	if(usetwopoints)
	{
		if(factofwholestep < lambda)
		{
			for(int i=0; i<3; i++)
			{
				posxyz[i] = Pos[i] + factofwholestep*dir1[i];
			}
		}
		else
		{
			double factminuslambd = factofwholestep-lambda;
			for(int i=0; i<3; i++)
			{
				posxyz[i] = xyz1[i] + factminuslambd*dir2[i];
			}
		}
	}
	else
	{
		for(int i=0; i<3; i++)
		{
			posxyz[i] = Pos[i] + factofwholestep*dir1[i];
		}
	}
	
	B->eval(t+factofwholestep*htry,posxyz,Bvec);
	return;	
}

double *Tracking::getLastPos()
{
	return Pos;
}

void Tracking::getInitialB(double t, double* Bvec)
{
	B->eval(t,Pos,Bvec);
	return;
}

void Tracking::initialize()
{
	double r = Rsquare+1.0;
	t = 0.0;
	while(r > Rsquare)
	{
		r = 0.0;
		for(int i=0; i<3; i++)
		{
			#pragma omp critical
			{
				Pos[i] = (2.*rand->uniform()-1.)*R;
			}
			r += Pos[i]*Pos[i];
		}
	}
}

void Tracking::stepDone()
{
	if(usetwopoints)
	{
		if(scaling >= 0.9999999)
		{
			for(int i=0; i<3; i++)
			{
				Pos[i] = xyz2[i];
			}
		}
		else if(scaling > lambda)
		{
			double sc = scaling - lambda;
			for(int i=0; i<3; i++)
			{
				Pos[i] = xyz1[i] + sc*dir2[i];
			}
		}
		else
		{
			for(int i=0; i<3; i++)
			{
				Pos[i] += scaling*dir1[i];
			}
		}
	}
	else
	{
		if(scaling >= 0.9999999)
		{
			for(int i=0; i<3; i++)
			{
				Pos[i] = xyz1[i];
			}
		}
		else
		{
			for(int i=0; i<3; i++)
			{
				Pos[i] += scaling*dir1[i];
			}
		}
	}
}

void Tracking::setScaling(double scal)
{
	scaling = scal;
}

void Tracking::makeTrack(double &time, double &ht)
{
	htry = ht;
	t = time;
	double sigma = sqrt2D*sqrt(htry);
	double r = 0.;
	scaling = 1.;
	for(int i=0; i<3; i++)
	{
		#pragma omp critical
		{
			dir1[i] = rand->gaussian(sigma);
		}
		xyz1[i] = Pos[i]+dir1[i];
		r += xyz1[i]*xyz1[i];
	}
	if(r > Rsquare)
	{
		usetwopoints = true;
		int NreflectionsOutside = -1;
		double a = 0.0;
		double b = 0.0;
		double rPosSquare = 0.0;
		for(int i=0; i<3; i++)
		{
			a += dir1[i]*dir1[i];
			b += Pos[i]*dir1[i];
			rPosSquare += Pos[i]*Pos[i];
		}
		if(a > Rsquare)
		{
			throw "Stepsize too big for diffusion!";
			exit(1);
		}
		lambda = (sqrt(b*b - a*(rPosSquare-Rsquare)) - b) / a;
		for(int i=0; i<3; i++)
		{
			xyz1[i] = Pos[i] + lambda*dir1[i];
		}
		r = Rsquare+1;
		while(r > Rsquare)
		{
			NreflectionsOutside++;
			r = 1.1;
			while(r > 1.)
			{
				r = 0.0;
				for(int i=0; i<3; i++)
				{
					#pragma omp critical
					{
						dir2[i] = 2*rand->uniform()-1.;
					}
					r += dir2[i]*dir2[i];
				}
			}
			r = sqrt(r);
			double sqrootofa = sqrt(a);
			for(int i=0; i<3; i++)
			{
				dir2[i] *= sqrootofa/r;
				xyz2[i] = xyz1[i] + (1-lambda)*dir2[i];
			}
			r = 0.0;
			for(int i=0; i<3; i++)
			{
				r += xyz2[i]*xyz2[i];
			}
		}
	}
	else
	{
		usetwopoints = false;
	}
}
