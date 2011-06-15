#include "tubegeometry.h"
#include "debug.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
using namespace std;

// Im Konstruktor wird ein File ausgelesen, das die Informationen über den Schlauch 
// enthält. Der Schlauch soll sich aus Kreisbogenstücken und geraden Stücken
// zusammensetzen. Der Konstruktor erzeugt aus dem File einen Vektor mit Segmenten
// des Schlauchs. Diese Segmente sind entweder Objekte der Klasse Csegment oder Lsegment.
Tubegeometry::Tubegeometry(Random *ran, std::string Tubefile, Parameters& theParameters)
	:Basegeometry(ran)
{
	ifstream tube;
	std::string dummy;
	tube.open(Tubefile.c_str());
	int numberofsegments;
	if(tube.is_open()) {
		string type;
		double s1,s2,s3,v1,v2,v3;
		tube >> numberofsegments;
		std::getline (tube, dummy);
		std::getline (tube, dummy);
// Liest Startvektor und Startachse der Schlauchachse aus.
		tube >> s1 >> s2 >> s3 >> v1 >> v2 >> v3;
		Threevector s (s1,s2,s3);
		Threevector v (v1,v2,v3);
		std::getline (tube, dummy);
		std::getline (tube, dummy);
		tube >> type;
// Das erste Segment wird erzeugt (Sonderbehandlung für erstes Segment).
// Wenn das 1. Segment ein Kreissegment ist wird noch der Winkel 't',
// der Radius 'r' des Bogens und der Normalenvektor 'n' des Kreises ausgelesen
// und ein Kreissegment Csegment erstellt.
		if (type == "c") {			
			double t,r,n1,n2,n3;
			tube >> t >> r >> n1 >> n2 >> n3;
			std::getline (tube, dummy);
			Threevector n(n1,n2,n3);	
			//Csegment C(s,v,n,r,t,theParameters);
			//Csegments.push_back(C);			
			//Segments.push_back(& Csegments.back());
			
			// TODO: memory leak
			Csegment *cs = new Csegment(s,v,n,r,t,theParameters);
			Segments.push_back(cs);
			std::cout << (*cs).toString() << std::endl;
		}
// Ist das erste Segment ein gerades Segment, wird die Länge 't' und nochmals die
// Richtung 'v' (die mit der zuerst ausgelesenen Startrichtung übereinstimmen sollte)
// ausgelesen und daraus ein gerades Segment Lsegment erzeugt.
		else if (type == "l") {
			double t,v1,v2,v3;
			tube >> t >> v1 >> v2 >> v3;
			std::getline (tube, dummy);
			Threevector v(v1,v2,v3);
			
			//Lsegment L(s,v,t,theParameters);
			//Lsegments.push_back(L);			
			//Segments.push_back(& Lsegments.back());
			
			// TODO: memory leak
			Lsegment *ls = new Lsegment(s,v,t,theParameters);
			Segments.push_back(ls);
			std::cout << (*ls).toString() << std::endl;
		}
		else{ 
			std::cerr << "File not good \n"; 			
		}

// Bei den anderen Einträgen wird ähnlich verfahren, nur dass der jeweilige Startvektor
// der Endvektor des vorherigen Segments ist und bei Kreissegmenten der Startachsenvektor
// auch der Endachsenvektor des vorherigen Segments ist.
		//while (tube.good()) {	
		for (int i = 0; i < numberofsegments-1; i++) {		
			tube >> type;
			if (type == "c") {
				Threevector s,v;
				double r, t, n1, n2, n3;
				s = Segments.back()->getposition(1);		
				v = Segments.back()->axis(1);
				tube >> t >> r >> n1 >> n2 >> n3;
				std::getline (tube, dummy);
				Threevector n(n1,n2,n3);					
				//Csegment C(s,v,n,r,t,theParameters);
				//Csegments.push_back(C);			
				//Segments.push_back(& Csegments.back());
			
				// TODO: memory leak
				Csegment *cs = new Csegment(s,v,n,r,t,theParameters);
				Segments.push_back(cs);
				std::cout << (*cs).toString() << std::endl;
			}
			else if (type == "l") {
								double 	t,v1,v2,v3;
				Threevector s;
				tube >> t >> v1 >> v2 >> v3;
				std::getline (tube, dummy);
				s = Segments.back()->getposition(1) ;
				Threevector v(v1,v2,v3);								
				//Lsegment L(s,v,t,theParameters);
				//Lsegments.push_back(L);			
				//Segments.push_back(& Lsegments.back());
			
				// TODO: memory leak
				Lsegment *ls = new Lsegment(s,v,t,theParameters);
				Segments.push_back(ls);
				std::cout << (*ls).toString() << std::endl;
			}
			else 
				std::cerr << "File not good \n";	
		}		
	}
	else 
		std::cerr << "Could not open file \n" ;
	tube.close(); 
	radiustube = theParameters.getDoubleParam("radiustube");
	std::cout << "Anzahl der Segmente: " << Segments.size()-1 << std::endl;
}

void Tubegeometry::initialize(Threevector &v, Threevector &x){
	Threevector standard;
	standard = Threevector(0,-1,0);
	debug << standard.toString() << std::endl;
	debug << Segments[0]->axis(0).toString() << std::endl;
	if (Segments[0]->axis(0).compare(standard)){
		Threevector a, b;
		double r, angle;
		a = Threevector (1,0,0);
		b = Threevector (0,0,1);
		r = fRandom->uniform() * radiustube; 
		angle = fRandom->uniform() * 360; 
		x = Segments[0]->startpoint() + r * cos(angle) * a + r * sin(angle) * b;
		v = Segments[0]->axis(0);
		debug << x.toString() << v.toString() << std::endl;
	}
	else {
		Threevector a, b, d;
		double r, angle;
		d = Segments[0]->axis(0);
		a = Threevector ((-1)*(d[1] + d[2])/d[0], 1, 1);
		b = a.normalized().cross(d.normalized());
		r = fRandom->uniform() * radiustube; 
		angle = fRandom->uniform() * 360; 
		x = Segments[0]->startpoint() + r * cos(angle) * a.normalized() + r * sin(angle) * b;
		v = Segments[0]->axis(0);
	}
}
// 'contains' prüft ob eine Position x im Schlauch enthalten ist oder nicht.
// Wenn ja, gibt sie den axialen Vektor der Projektion der Position auf die Schlauchachse
// zurück. Wenn nicht, dann gibt sie als axialen Vektor den Nullvektor zurück.

//Threevector Tubegeometry::contains(const Threevector &x, double &rootstart)  {
Threevector Tubegeometry::contains(const Threevector &x)  {
	Threevector v(0,0,0);
	for (size_t i = 0; i < Segments.size(); i++) {	
		//double rootstartsegment;
		bool inorout;
		Threevector axis;
		//rootstartsegment = rootstart - (double)i;
		//axis = Segments[i]->segmentcontains(x, rootstartsegment);
		axis = Segments[i]->segmentcontains(x);
		inorout = axis.compare(v);
		if (inorout == false){
			//rootstart = rootstartsegment + (double)i;
			return axis;
		}
	}
	return v;
}

bool Tubegeometry::lastsegmentcontains(const Threevector &x) {
	//double rootstartsegment;
	Threevector v(0,0,0);
	Segment* lastsegment = Segments.back();
	//rootstartsegment = 0;
	//if (lastsegment->segmentcontains(x,rootstartsegment).compare(v))
	if (lastsegment->segmentcontains(x).compare(v))
		return false;
	else 
		return true;
}
