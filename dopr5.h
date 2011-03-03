#ifndef DOPR5_H
#define DOPR5_H

#include "derivatives.h"
#include "tracking.h"

class Dopr5
{
	public:
		Dopr5(double, unsigned int, double*, double*, double, Derivatives *, double, double, bool);
		~Dopr5();
		void reset(double, double*, double*, double);
		void step();
		void dy(const double);
		double error();
		bool success(const double, double&);
		void prepare_dense(const double);
		double dense_out(const int, const double, const double);
		
		inline double getHdid(){return hdid;};
		inline double getT(){return t;};
	
	private:	
		double told, hdid, htry, hnext;
		const unsigned int n;
		double *y, *yout, *dydx, *yerr;
		double t;
		double *k2,*k3,*k4,*k5,*k6;
		double *rcont1, *rcont2, *rcont3, *rcont4, *rcont5;
		double *dydxnew;
		Derivatives *derivatives;
		double atol, rtol, EPS, errold;
		Tracking *tracker;
		bool reject;
		bool dense;
};
#endif
