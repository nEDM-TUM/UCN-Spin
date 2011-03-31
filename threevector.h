#ifndef THREEVECTOR_H
#define THREEVECTOR_H

#include <vector>

class Threevector
{
	public:
		Threevector(double x, double y, double z);
		Threevector();
		Threevector(double *v);
		void setX(double x){fVec[0] = x;};
		void setY(double y){fVec[1] = y;};
		void setZ(double z){fVec[2] = z;};
		double getX(){return fVec[0];};
		double getY(){return fVec[1];};
		double getZ(){return fVec[2];};
		double &operator[](const int i) {return fVec[i];};
		Threevector &operator=(const Threevector &v) {
			for (int i = 0; i < 3; i++)
				fVec[i] = v.fVec[i];
			return *this;
		};
		Threevector &operator+=(const Threevector &v) {
			for (int i = 0; i < 3; i++)
				fVec[i] += v.fVec[i];
			return *this;
		}
		friend Threevector operator+(const Threevector &left, const Threevector &right);
		friend double operator*(const Threevector &left, const Threevector &right);
		friend Threevector operator*(const double &left, const Threevector &right);
		friend Threevector operator*(const Threevector &left, const double &right) {
			return right * left;
		};
		friend bool operator==(const Threevector &left, const Threevector &right);
		friend bool operator!=(const Threevector &left, const Threevector &right);
		double mag();
		double magsquare();

	protected:
		double fVec[3];
};

#endif
