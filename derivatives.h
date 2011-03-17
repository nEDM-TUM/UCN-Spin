#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "tracking.h"

class Derivatives
{
	public:
		Derivatives(Tracking*);
		void operator()(const double, const double[], double[]);
		void eval(const double, const double[], double[]);
		void evalInitial(const double, const double[], double[]);
		inline Tracking *getTracker(){return tracker;};
	private:
		Tracking *tracker;
};
#endif
