#include "dipolefield.h"
#include "globals.h"
#include <iostream>
#include <math.h>

using namespace std;

Dipolefield::Dipolefield(double B0)
: factor(-3.9293341e-34), omegalarmor(-PRECES_N*B0)
{
	
}

void Dipolefield::getField(double t, double *Bout, double *r, double *P)
{
	double co = cos(omegalarmor*t);
	double si = sin(omegalarmor*t);
	double Prot[3] = {P[0]*co - P[1]*si, P[0]*si + P[1]*co, P[2]};
	double dist[3] = {r[0],r[1],r[2]};
	double dist_mag = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
	for(int i=0; i<3; i++)
	{
		dist[i] /= dist_mag;
	}
	double DistFact = factor/(dist_mag*dist_mag*dist_mag);
	double ScalProd = 3*(Prot[0]*dist[0]+Prot[1]*dist[1]+Prot[2]*dist[2]);
	Bout[0] = DistFact*(ScalProd*dist[0]-Prot[0]);
	Bout[1] = DistFact*(ScalProd*dist[1]-Prot[1]);
	Bout[2] = DistFact*(ScalProd*dist[2]-Prot[2]);
}
