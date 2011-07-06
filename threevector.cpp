#include "debug.h"
#include "threevector.h"
#include <cmath>
#include <sstream>

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

Threevector::Threevector(const Threevector &v)
{
	//debug << "Copy constructor called!" << std::endl;
	for(int i=0; i<3; i++)
		fVec[i] = v.fVec[i];
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

bool Threevector::compare(const Threevector &x){
	Threevector y = *this;
	for (int i = 0; i < 3; i++){
		if (fabs(y[i]-x[i]) > 1e-10)
			return false;
	}
	return true;
		
}

double Threevector::mag() const
{
	double a=0.0;
	for(int i=0; i<3; i++)
	{
		a += fVec[i]*fVec[i];
	}
	return sqrt(a);
}

double Threevector::magsquare() const
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

Threevector Threevector::normalized()
{
	Threevector temp(*this);
	temp.normalize();
	return temp;
}

std::string Threevector::toString() const {
	std::ostringstream o;
	o.precision(8);
	o << "(" << fVec[0] << " " << fVec[1] << " " << fVec[2] << ")";
	return o.str();
}

Threevector Threevector::cross (const Threevector &x) const {
	Threevector c;
	c[0]=fVec[1]*x[2]-fVec[2]*x[1];
	c[1]=fVec[2]*x[0]-fVec[0]*x[2];
	c[2]=fVec[0]*x[1]-fVec[1]*x[0];
	return c;
}
