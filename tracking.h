#ifndef TRACKING2_H
#define TRACKING2_H

#include "random.h"
#include "bfield.h"

class Tracking
{
	public:
		Tracking(Random*, Bfield*, double, double);
		~Tracking();
		void initialize();
		void getB(double, double*);
		double *getLastPos();
		void getInitialB(double, double*);
		void makeTrack(double&,double&);
		void setScaling(double);
		void stepDone();
	private:
		Bfield *B;
		Random *rand;
		double *Pos, *xyz1, *xyz2, *dir1, *dir2, *posxyz;
		double sqrt2D, lambda, Rsquare, R, scaling, t;
		bool usetwopoints;
		double htry;
};
#endif
