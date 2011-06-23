#ifndef THREEVECTOR_H
#define THREEVECTOR_H

#include "debug.h"
#include <vector>
#include <string>
#include <iostream>

class Threevector
{
	public:
		Threevector(double x, double y, double z);
		Threevector();
		Threevector(double *v);
		Threevector(const Threevector &t);
		void setX(double x){fVec[0] = x;};
		void setY(double y){fVec[1] = y;};
		void setZ(double z){fVec[2] = z;};
		double getX(){return fVec[0];};
		double getY(){return fVec[1];};
		double getZ(){return fVec[2];};
		double &operator[](const int i) {return fVec[i];};
		double operator[](const int i) const {return fVec[i];};
		Threevector &operator*=(const double a) {
			for (int i = 0; i < 3; i++)
				fVec[i] *= a;
			return *this;
		};
		Threevector &operator=(const Threevector &v) {
			for (int i = 0; i < 3; i++)
				fVec[i] = v.fVec[i];
			return *this;
		};
		Threevector &operator+=(const Threevector &v) {
			for (int i = 0; i < 3; i++)
				fVec[i] += v.fVec[i];
			return *this;
		};
		Threevector operator-() const {
			Threevector ret;
			for (int i = 0; i < 3; i++)
				ret[i] = -fVec[i];
			return ret;
		};

		friend Threevector operator+(const Threevector &left, const Threevector &right);
		friend Threevector operator-(const Threevector &left, const Threevector &right);
		friend double operator*(const Threevector &left, const Threevector &right);
		friend Threevector operator*(const double &left, const Threevector &right);
		friend Threevector operator*(const Threevector &left, const double &right) {
			return right * left;
		};
		friend bool operator==(const Threevector &left, const Threevector &right);
		friend bool operator!=(const Threevector &left, const Threevector &right);
		double mag() const;
		double magsquare() const;
		void normalize();
		std::string toString() const;
		Threevector cross (const Threevector &x) const ;
	protected:
		double fVec[3];
};

std::ostream& operator<<(std::ostream &o, const Threevector &v);
std::istream& operator>>(std::istream &i, Threevector &v);

#endif
