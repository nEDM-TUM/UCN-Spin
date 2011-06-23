#ifndef _BFIELD_H
#define _BFIELD_H

#include "parameters.h"
#include "basetracking.h"
#include "threevector.h"

class Bfield
{
	public:
		Bfield(Basetracking *tracker) : fTracker(tracker) {};
		virtual ~Bfield() {};
		Threevector operator()(const double time) const { return eval(time); };
		virtual Threevector eval(const double time) const = 0;

	protected:
		Basetracking *fTracker;
};

#endif // _BFIELD_H
