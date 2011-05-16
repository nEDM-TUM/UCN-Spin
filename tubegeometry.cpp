#include "tubegeometry.h"
#include <fstream>
#include <iostream>
#include <vector>

Tubegeometry::Tubegeometry(Random *ran, std:string Tubefile){
	ifstream tube;
	tube.open(Tubefile.c_str());
	vector<Segment> segments;
	if(tube.is_open()) {
		string type;
		double s1,s2,s3,v1,v2,v3;
		getline (tube, comment1);
// Liest Startvektor und Startachse der Schlauchachse aus.
		tube >> s1 >> s2 >> s3 >> v1 >> v2 >> v3;
		Threevector s (s1,s2,s3);
		Threevector v (v1,v2,v3);
		getline (tube, comment2);
		tube >> type;
// Das erste Segment wird erzeugt (Sonderbehandlung für erstes Segment).
// Wenn das 1. Segment ein Kreissegment ist wird noch der Winkel 't',
// der Radius 'r' des Bogens und der Normalenvektor 'n' des Kreises ausgelesen
// und ein Kreissegment Csegment erstellt.
		if (type == "c") {				
			double t,n1,n2,n3;
			tube >> t >> r >> n1 >> n2 >> n3;
			Threevector n(n1,n2,n3);	
			Csegment O (s,v,n,r,t);
			segments.push_back(O); 		
		}
// Ist das erste Segment ein gerades Segment, wird die Länge 't' und nochmals die
// Richtung 'v' (die mit der zuerst ausgelesenen Startrichtung übereinstimmen sollte)
// ausgelesen und daraus ein gerades Segment Lsegment erzeugt.
		else if (type == "l") {
			double t,v1,v2,v3;
			tube >> t >> v1 >> v2 >> v3;
			Threevector v(v1,v2,v3);
			Lsegment O(s,v,t);
			segments.push_back(O);			
		}
		else 
			cerr << "File not good \n"; 			

// Bei den anderen Einträgen wird ähnlich verfahren, nur dass der jeweilige Startvektor
// der Endvektor des vorherigen Segments ist und bei Kreissegmenten der Startachsenvektor
// auch der Endachsenvektor des vorherigen Segments ist.
		while (tube.good()) {			
			tube >> type;
			if (type == "c") {
				Threevector s,v;
				double r, t, n1, n2, n3;
				s = segments.back().getposition(1);		
				v = segments.back().axis(1);
				tube >> t >> r >> n1 >> n2 >> n3;
				Threevector n(n1,n2,n3);					
				Csegment O (s,v,n,r,t);
				segments.push_back(O);
			}
			else if (type == "l") {
				double 	t,v1,v2,v3;
				tube >> t >> v1 >> v2 >> v3;
				s = segments.back().getposition(1) ;
				Threevector v(v1,v2,v3);								
				Lsegment O(s,v,t);
				segments.push_back(O);
			}
			else 
				cerr << "File not good \n";			
			}
	else 
		cerr << "Could not open file \n" ;
	tube.close(); 
	radiustube = theParameters.getDoubleParam("radiustube");
}

void TubeGeometry::initialize(Threevector &v, Threevector &x){
	standard = Threevector(0,-1,0);
	if (segments[0].axis(0) == standard){
		Threevector a, b;
		double r, angle;
		a = Threevector (1,0,0);
		b = Threevector (0,0,1);
		r = random->uniform() * radiustube; 
		angle = random->uniform() * 360; 
		x = segments[0].start() + r * cos(angle) * a + r sin(angle) * b;
		v = segments[0].axis(0);
	}
	else {
		Threevector a, b, d;
		double r, angle;
		d = segments[0].axis(0);
		a = Threevector ((-1)*(d[1] + d[2])/d[0], 1, 1);
		b = a.normalize().cross(d.normalize());
		r = random->uniform() * radiustube; 
		angle = random->uniform() * 360; 
		x = segments[0].start() + r * cos(angle) * a.normalize() + r sin(angle) * b;
		v = segments[0].axis(0);
	}
}


Threevector TubeGeometry::contains(const Threevector &x, double &rootstart)  {
	Threevector v;
	v = Threevector();
// Testet in jedem Segment, ob der Punkt enthalten ist oder nicht. Sobald ein Segment
// einen Axialenvektor ungleich dem Nullvektor v zurückgibt, wird dieser ausgegeben.
// Gleichzeitig wird dann der Wert 'rootstart' verändert. Er ist dann die momentane
// Projektion der Position auf die Schlauchachse und soll im nächsten Schritt als
// Startwert für die Nullstellensuche verwendet werden.  
	for (vector<Segments>::const_iterator it = segments.begin(); it !=segments_end(); it++) {
		double rootstartsegment;
		rootstartsegment = rootstart - it;
		if (*it.segmentcontains(x, trootstartsegment) != v){
			rootstart = rootstartsegment + it
			return *it.segmentcontains(x, rootstartsegment);
		}
	}
	return v;
}

bool TubeGeometry::lastsegmentcontains(const Threevector &x) {
	Segment lastsegment;
	double rootstartsegment;
	lastsegment = segment.back();
	rootstartsegment = 0
	if (lastsegment.segmentcontains(x,rootstartsegment) == true)
		return true;
	else 
		return false;
}
