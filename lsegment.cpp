# include "threevector.h"
# include "roots.h"
# include <cmath>
# include "lsegment.h"


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
Threevector Lsegment::segmentcontains(const Threevector &x, double &rootstartsegment) {
	double projection;
	Threevector radial;
	Threevector axial;
	Threevector v(0,0,0);
	projection = (x + (-1 * start)) * direction;
	if (projection < 0)
		return v;
	else if (projection > t_max)
		return v;
	else {
		radial = x + (-1) * (projection * direction);
		if (radial.magsquare() > radiustube * radiustube) 
			return v; 
		else {
			return axis();
			rootstartsegment = projection / t_max;
		}
	}
}  
