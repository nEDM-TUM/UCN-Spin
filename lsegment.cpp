#include "threevector.h"
#include  "debug.h"
#include "roots.h"
#include <cmath>
#include "lsegment.h"
#include <sstream>
#include <string>

Lsegment::Lsegment(Threevector s, Threevector v, double t, Parameters& theParameters) :
	start(s), direction(v.normalized()), t_max(t), radiustube(theParameters.getDoubleParam("radiustube"))
{
}

Threevector Lsegment::getposition (double tau) {
	Threevector pos;
	pos = start + tau * t_max * direction;
	return pos;
}

Threevector Lsegment::axis(double tau) {
	return axis();
}

Threevector Lsegment::axis () {
	return direction;
}

Threevector Lsegment::startpoint() {
	return start;
}

// 'rootstartsegment' wird während der Funktion verändert.
//Threevector Lsegment::segmentcontains(const Threevector &x, double &rootstartsegment) {
Threevector Lsegment::segmentcontains(const Threevector &x) {
	double projection;
	Lsegment L = *this;
	Threevector radial, axial;
	Threevector v(0,0,0);
	projection = (x + (-1 * L.startpoint())) * direction;
	
	if (projection < 0){
		return v;
	}
	else if (projection > t_max){
		return v;
	}
	else {
		radial = x + (-1) * L.startpoint() + (-1) * (projection * direction);
		if (radial.magsquare() > radiustube * radiustube) {
			return v; 
		}
		else {
			//rootstartsegment = projection / t_max;
			return L.axis().normalized();
		}
	}
}  

std::string Lsegment::toString() const {
	std::ostringstream o;
	o << "L-Segment: " << "Startpunkt = " << (*this).start.toString()  << ";		Richtung = " << (*this).direction.toString() << ";	Länge = " << t_max;
	return o.str();
}
