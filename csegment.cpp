# include "threevector.h"
# include "roots.h"
# include <cmath>
# include "csegment.h"
# include <sstream>
# include <string>

Csegment::Csegment(Threevector s, Threevector v, Threevector n, double r, double t, Parameters& theParameters) :
	radius(r), radiustube(theParameters.getDoubleParam("radiustube")), t_max(t), start(s), normal(n), b(v.normalized()), a(b.cross(n.normalized())),
	centre(start + (-radius)*a)
{
	position = Threevector (0.0, 0.0, 0.0);
}

Threevector Csegment::getposition (double tau) { 
	Threevector pos;
	pos = centre + radius*cos(tau * t_max)*a + radius*sin(tau * t_max)*b;
	return pos;
}
 
Threevector Csegment::axis(double tau) {
	Threevector ax;
	ax = (-1)*sin(tau * t_max)*a + cos(tau * t_max)*b;
	return ax;
}


Threevector Csegment::segmentcontains(const Threevector &x) {
	double mindistancesquare, tau_0;
	Threevector c, v(0,0,0);
	position = x; 
	if (derivdist(0.0)*derivdist(1.0) > 0.0){
		return v;
	}
	tau_0 = Roots::safeNewton(this,&Csegment::derivdist, &Csegment::secderivdist, 0.0, 1.0, 0.0000001);
	if (secderivdist(tau_0) < 0){
		return v;
	}
	else {
		c = x + (-1) * getposition(tau_0);
		mindistancesquare = c.magsquare();
		if (mindistancesquare > radiustube*radiustube){
			return v;
		}
		else {
			return axis(tau_0);
		}
	}
}


double Csegment::derivdist(double tau) {
	double derivdistance = 0;
	Threevector k, l;
	k = Threevector();
	l = Threevector();
	for (int i = 0; i < 3; i++) {
		k[i] = radius * (cos(tau * t_max)*a[i] + sin(tau * t_max)*b[i]);
		l[i] = radius * (-sin(tau * t_max)*a[i] + cos(tau * t_max)*b[i]);
	}
	for (int i = 0; i<3; i++) {
		derivdistance = derivdistance + 2*(centre[i] + k[i] - position[i]) * l[i] * t_max; 
	}
	return derivdistance;
}

double Csegment::secderivdist(double tau ) {
	double secderivdistance = 0;
	Threevector k, l;
	k = Threevector();
	l = Threevector();
	for (int i = 0; i<3; i++) {
		k[i] = radius * (cos(tau * t_max)*a[i] + sin(tau * t_max)*b[i]);
		l[i] = radius * (-sin(tau * t_max)*a[i] + cos(tau * t_max)*b[i]);
	}
	for (int i  =0; i<3; i++) 
		secderivdistance = secderivdistance + 2*(centre[i] + k[i] - position[i]) * (-k[i]) * t_max * t_max + 2*(l[i]*l[i]) * t_max * t_max; 
	return secderivdistance;
}

std::string Csegment::toString() const {
	std::ostringstream o;
	o << "C-Segment: " << "r(tau) = " << centre.toString() << " + " << radius << " * " << "cos(tau * " << t_max << ") * " << a.toString() << " + " << radius << " * " << "sin(tau * " << t_max << ")* " << b.toString();
	return o.str();
}

std::string Csegment::toMathematica() const{
	std::ostringstream o;
	o << "{" << centre[0] << " + " << radius << " * Cos[x * " << t_max << "] * " << a[0] << " + " << radius << " * Sin[x * " << t_max << "] * " << b[0] << ", " << centre[1] << " + " << radius << " * Cos[x * " << t_max << "] * " << a[1] << " + " << radius << " * Sin[x * " << t_max << "] * " << b[1] << ", " << centre[2] << " + " << radius << " * Cos[x * " << t_max << "] * " << a[2] << " + " << radius << " * Sin[x * " << t_max << "] * " << b[2] << "}"; 
	return o.str();
}
