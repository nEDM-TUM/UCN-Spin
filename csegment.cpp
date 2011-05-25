# include "threevector.h"
# include "roots.h"
# include <cmath>
# include "csegment.h"

Csegment::Csegment(Threevector s, Threevector v, Threevector n, double r, double t, Parameters& theParameters) :
	radius(r), radiustube(theParameters.getDoubleParam("radiustube")), t_max(t), start(s), b(v.normalized()), a(b.cross(n.normalized())),
	centre(start + (-radius)*a)
{
//	radius = r;
//	t_max = t;
//	start = s;
//	b = v.normalize();
//	a = b.cross(n.normalize());
//	centre = start + (-radius)*a;
//	radiustube = theParameters.getDoubleParam("radiustube");
}

Threevector Csegment::getposition (double tau) { 
	Threevector pos;
	pos = centre + radius*cos(tau * t_max)*a + radius*sin(tau * t_max)*b;
	return pos;
}

// Gibt die Achsenrichtung des Schlauchs an einem Punkt auf der Achse wieder. 
Threevector Csegment::axis(double tau) {
	Threevector ax;
	ax = (-1)*sin(tau * t_max)*a + cos(tau * t_max)*b;
	return ax;
}

Threevector Csegment::startpoint() {
	return start;
}

// Überprüft ob ein Vektor x im Schlauch enthalten ist. Dazu muss mit rootNewton ein
// Minimum gesucht werden. Der Startwert 'rootstartsegment' der Nullstellensuche
// wird während der Funktion verändert.
Threevector Csegment::segmentcontains(const Threevector &x, double &rootstartsegment) {
	double mindistancesquare, tau_0;
	Threevector c, v(0,0,0);
	tau_0 = rootNewton(x, rootstartsegment);
	if (tau_0 == 2)
		return v;
	else if (secderivdist(x,tau_0) < 0)
		return v;
	else {
		c = x + (-1) * getposition(tau_0);
		mindistancesquare = c.magsquare();
		if (mindistancesquare > radiustube*radiustube)
				return v;
		else {
			rootstartsegment = tau_0;
			return axis(tau_0);
		}
	}
}

double Csegment::rootNewton(const Threevector &x, double rootstartsegment){
	double root; 
	if (rootstartsegment < 0)
		root = 0;
	else if (rootstartsegment > 1)
		root = 1;
	else 
		root = rootstartsegment; 
	if (derivdist(x,0)*derivdist(x,1) > 0){
		return 2;
	}
	else {
		while (fabs(derivdist(x,root)) > 0.000001) {
			root = root - (derivdist(x,root) / secderivdist(x,root));
		}
	}
	return root;
}

double Csegment::derivdist(const Threevector &x, double tau) {
	double derivdistance = 0;
	Threevector k, l;
	k = Threevector();
	l = Threevector();
	for (int i = 0; i < 3; i++) {
		k[i] = radius * (cos(tau * t_max)*a[i] + sin(tau * t_max)*b[i]);
		l[i] = radius * (-sin(tau * t_max)*a[i] + cos(tau * t_max)*b[i]);
	}
	for (int i=0; i<3; i++) 
		derivdistance = derivdistance + 2*(centre[i] + k[i] - x[i]) * l[i]; 
	return derivdistance;
}

double Csegment::secderivdist(const Threevector &x, double tau ) {
	double secderivdistance;
	Threevector k, l;
	k = Threevector();
	l = Threevector();
	for (int i = 0; i<3; i++) {
		k[i] = radius * (cos(tau * t_max)*a[i] + sin(tau * t_max)*b[i]);
		l[i] = radius * (-sin(tau * t_max)*a[i] + cos(tau * t_max)*b[i]);
	}
	secderivdistance = 0;
	for (int i=0; i<3; i++) 
		secderivdistance = secderivdistance + 2*(centre[i] + k[i] - x[i]) * (-k[i]) + 2*(l[i]*l[i]); 
	return secderivdistance;
}
