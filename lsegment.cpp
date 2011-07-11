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
	return direction;
}


Threevector Lsegment::segmentcontains(const Threevector &x) {
	double projection;
	Threevector radial, axial;
	Threevector v(0,0,0);
	projection = (x + (-1 * start)) * direction;
	
	if (projection < 0){
		return v;
	}
	else if (projection > t_max){
		return v;
	}
	else {
		radial = x + (-1) * start + (-1) * (projection * direction);
		if (radial.magsquare() > radiustube * radiustube) {
			return v; 
		}
		else {
			return direction;
		}
	}
}  

std::string Lsegment::toString() const {
	std::ostringstream o;
	o << "L-Segment: " << "r(tau) = " << start.toString() << " + " << "tau"  << " * "  << t_max  << " * " << direction.toString();
	return o.str();
}

std::string Lsegment::toMathematica() const {
	std::ostringstream o;
	o << "{" << start[0] << " + " << "x * " << t_max << " * " << direction[0] << ", " << start[1] << " + " << "x * " << t_max << " * " << direction[1] << ", "  << start[2] << " + " << "x * " << t_max << " * " << direction[2] << "}";
	return o.str();
}
