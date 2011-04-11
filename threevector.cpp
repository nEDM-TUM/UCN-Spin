#include "threevector.h"
#include <cmath>

Threevector::Threevector(double x, double y, double z)
{
	fVec[0] = x;
	fVec[1] = y;
	fVec[2] = z;
}

Threevector::Threevector()
{
	for(int i=0; i<3; i++)
	{
		fVec[i] = 0.0;
	}
}

Threevector::Threevector(double *v)
{
	for(int i=0; i<3; i++)
	{
		fVec[i] = v[i];
	}
}

Threevector operator+(const Threevector &left, const Threevector &right)
{
	Threevector a;
	for(int i=0; i<3; i++)
	{
		a.fVec[i] = left.fVec[i] + right.fVec[i];
	}
	return a;
}

double operator*(const Threevector &left, const Threevector &right)
{
	double a=0.0;
	for(int i=0; i<3; i++)
	{
		a += left.fVec[i]*right.fVec[i];
	}
	return a;
}

Threevector operator*(const double &left, const Threevector & right)
{
	Threevector a;
	for(int i=0; i<3; i++)
	{
		a.fVec[i] = left*right.fVec[i];
	}
	return a;
}


bool operator==(const Threevector &left, const Threevector &right)
{
	bool a = true;
	for(int i=0; i<3; i++)
	{
		a = a && (left.fVec[i] && right.fVec[i]);
	}
	return a;
}

bool operator!=(const Threevector &left, const Threevector &right)
{
	return !(left==right);
}

double Threevector::mag()
{
	double a=0.0;
	for(int i=0; i<3; i++)
	{
		a += fVec[i]*fVec[i];
	}
	return sqrt(a);
}

double Threevector::magsquare()
{
	double a=0;
	for(int i=0; i<3; i++)
	{
		a += fVec[i]*fVec[i];
	}
	return a;
}

void Threevector::normalize()
{
	const double div = mag();

	for(int i=0; i<3; i++)
		fVec[i] /= div;
}
