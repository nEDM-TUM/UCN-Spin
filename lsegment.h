#ifndef LSEGMENT_H
#define LSEGMENT_H
#include "segment.h"
#include "threevector.h"

class Lsegment : public Segment {
	public:
		Lsegment (Threevector s, Threevector v, double t, Parameters& theParameters);
		Threevector getposition (double tau);
		Threevector axis(double tau);
		Threevector startpoint();
// Achtung 'rootstartsegmet' wird her verändert!!!
		//Threevector segmentcontains(const Threevector &x, double &rootstartsegment);
		Threevector segmentcontains(const Threevector &x);
		std::string toString() const;
	
	private:
		Threevector	start;
		Threevector direction;
		double t_max, radiustube;
		Threevector axis();

};

#endif
