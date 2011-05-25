#ifndef SEGMENT_H
#define SEGMENT_H
#include "threevector.h"
#include "parameters.h"

class Segment {
	public:
		virtual Threevector getposition (double tau)=0;
		virtual Threevector segmentcontains(const Threevector &x, double &rootstartsegment){return Threevector();};
		virtual Threevector axis(double tau){return Threevector();};
		virtual Threevector startpoint(){return Threevector();};
};

#endif
