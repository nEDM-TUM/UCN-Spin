# include "threevector.h"
# include "roots.h"
# include <cmath>

Csegment::Csegment(Threevector s, Threevector v, Threevector n, double r, double t) {
	radius = r;
	t_max = t;
	start = s
	b = v.normalize();
	a = b.cross(n.normalize);
	centre = start + (-radius)*a;
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
	double mindistancesquare;
	double tau_0;
	Threevector c;
	Threevector v;
	v = Threevector();
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
		while (abs(derivdist(x,root)) > 0.00001) {
			root = root - (derivdist(x,root) / secderivdist(x,root));
		}
	}
	return root;
}

double Csegment::derivdist(const Threevector &x, double tau) {
	double derivdistance;
	double k
	double l
	k = radius * (cos(tau * t_max)*a[i] + sin(tau * t_max)*b[i]);
	l = radius * (-sin(tau * t_max)*a[i] + cos(tau * t_max)*b[i]);
	derivdistance = 0
	for (int i=0; i<3; i++) {
		derivdistance = derivdistance + 2*(centre[i] + k - x[i]) * l; 
	return derivdistance;
}

double Csegment::secderivdist(const Threevector &x, double t ) {
	double secderivdistance;
	double k
	double l
	k = radius * (cos(tau * t_max)*a[i] + sin(tau * t_max)*b[i]);
	l = radius * (-sin(tau * t_max)*a[i] + cos(tau * t_max)*b[i]);
	secderivdistance = 0
	for (int i=0; i<3; i++) {
		secderivdistance = secderivdistance + 2*(centre[i] + k - x[i]) * (-k) + 2*(l*l); 
	return secderivdistance;
}

Lsegment::Lsegment(Threevector s, Threevector v, double t) {
	start = s;
	direction = v.normalize();
	t_max = t;
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
	Threevector v;
	v = Threevector();
	projection = (x + (-1 * start)) * direction;
	if (projection < 0)
		return v;
	else if (projection > t_max)
		return v;
	else {
		radial = x + (-1) * (projection * direction);
		if (rad.magsquare() > radiustube * radiustube) 
			return v; 
		else 
			return axis();
			rootstartsegment = projection / t_max
		}  
}

