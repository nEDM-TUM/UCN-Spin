#ifndef CSEGMENT_H
#define CSEGMENT_H
#include "segment.h"
#include "threevector.h"

class Csegment : public Segment {
	public:
		Csegment (Threevector s, Threevector v, Threevector n, double r, double t, Parameters& theParameters);
		Threevector getposition (double tau);
		Threevector axis(double tau);
		Threevector startpoint();
// Achtung 'rootstartsegment' wird hier verändert!!
		//Threevector segmentcontains(const Threevector &x, double &rootstartsegment);
		Threevector segmentcontains(const Threevector &x);
		std::string toString() const;
	
	private:
		double derivdist( double tau);
		double secderivdist( double tau);
		double rootNewton(const Threevector &x, double rootstartsegment);
		double radius, radiustube, t_max;
		Threevector start, b, a, centre;
		Threevector position; 
};

#endif
