#include "tubegeometry.h"
#include "debug.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
using namespace std;


Tubegeometry::Tubegeometry(Random *ran, std::string Tubefile, std::string Tubefilemathematica,  Parameters& theParameters)
	:Basegeometry(ran)
{
	ifstream tube; 
	fstream tubemathematica;
	std::string dummy;
	tube.open(Tubefile.c_str());
	tubemathematica.open(Tubefilemathematica.c_str());
	int numberofsegments;
	
	if(tube.is_open() && tubemathematica.is_open()) {
		string type;
		double s1,s2,s3,v1,v2,v3;
		std::getline (tube,dummy);
		tube >> numberofsegments;
		std::getline (tube, dummy);
		std::getline (tube, dummy);
		tube >> s1 >> s2 >> s3 >> v1 >> v2 >> v3;
		Threevector s (s1,s2,s3);
		Threevector v (v1,v2,v3);
		std::getline (tube, dummy);
		std::getline (tube, dummy);
		tube >> type;
	
		if (type == "c") {			
			double t,r,n1,n2,n3;
			tube >> t >> r >> n1 >> n2 >> n3;
			std::getline (tube, dummy);
			Threevector n(n1,n2,n3);	
			
			// TODO: memory leak
			Csegment *cs = new Csegment(s,v,n,r,t,theParameters);
			Segments.push_back(cs);
				#pragma omp master
                {
					std::cout << (*cs).toString() << std::endl;
                }
				tubemathematica << "Table" << Segments.size() << " = Table[" << (*cs).toMathematica() << ", {x, 0., 1, 0.1}]" << std::endl;
		}
		
		else if (type == "l") {
			double t,v1,v2,v3;
			tube >> t >> v1 >> v2 >> v3;
			std::getline (tube, dummy);
			Threevector v(v1,v2,v3);
			
			// TODO: memory leak
			Lsegment *ls = new Lsegment(s,v,t,theParameters);
			Segments.push_back(ls);
				#pragma omp master
                {
					std::cout << (*ls).toString() << std::endl;
                }
            tubemathematica << "Table" << Segments.size() << " = Table[" << (*ls).toMathematica() << ", {x, 0., 1, 0.1}]" << std::endl;
		}
		else{ 
			std::cerr << "File not good \n"; 			
		}

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
			
				// TODO: memory leak
				Csegment *cs = new Csegment(s,v,n,r,t,theParameters);
				Segments.push_back(cs);
					#pragma omp master
					{
						std::cout << (*cs).toString() << std::endl;
                    }
                 tubemathematica << "Table" << Segments.size() << " = Table[" << (*cs).toMathematica() << ", {x, 0., 1, 0.1}]" << std::endl;
			}
			
			else if (type == "l") {
				double 	t,v1,v2,v3;
				Threevector s;
				tube >> t >> v1 >> v2 >> v3;
				std::getline (tube, dummy);
				s = Segments.back()->getposition(1) ;
				Threevector v(v1,v2,v3);								
			
				// TODO: memory leak
				Lsegment *ls = new Lsegment(s,v,t,theParameters);
				Segments.push_back(ls);
					#pragma omp master
                    {
						std::cout << (*ls).toString() << std::endl;
                    }
                tubemathematica << "Table" << Segments.size() << " = Table[" << (*ls).toMathematica() << ", {x, 0., 1, 0.1}]" << std::endl;
			}
			else 
				std::cerr << "File not good \n";	
		}		
	}
	else 
		std::cerr << "Could not open file \n" ;
	tube.close();
	tubemathematica.close(); 
	radiustube = theParameters.getDoubleParam("radiustube");
	 #pragma omp master
	{                        
		std::cout << "Anzahl der Segmente (inkl Schlusssegment): " << Segments.size() << std::endl;
	}
}

void Tubegeometry::initialize(Threevector &v, Threevector &x){
		Threevector a, b, c, d, e;
		double r, angle;
		e = Threevector (0., 0., 0.);
		d = (Segments[0]->axis(0)).normalized();
		c = Threevector (0., 0., 1.);
		a = c.cross(d);
		if (a.compare(e) == true){
			c = Threevector (1., 0., 0.);
			a = c.cross(d);
			r = fRandom->uniform() * radiustube; 
			angle = fRandom->uniform() * 360; 
			x = Segments[0]->getposition(0) + r * cos(angle) * a.normalized() + r * sin(angle) * b;
		}			
		else {
			b = a.cross(e);
			r = fRandom->uniform() * radiustube; 
			angle = fRandom->uniform() * 360; 
			x = Segments[0]->getposition(0) + r * cos(angle) * a.normalized() + r * sin(angle) * b;
		}
	v = Segments[0]->axis(0);	
}


Threevector Tubegeometry::contains(const Threevector &x)  {
	Threevector v(0,0,0);
	for (size_t i = 0; i < Segments.size(); i++) {	
		bool inorout;
		Threevector axis;
		axis = Segments[i]->segmentcontains(x);
		inorout = axis.compare(v);
		if (inorout == false){
			return axis;
		}
	}
	return v;
}

bool Tubegeometry::lastsegmentcontains(const Threevector &x) {
	Threevector v(0,0,0);
	Segment* lastsegment = Segments.back();
	if (lastsegment->segmentcontains(x).compare(v))
		return false;
	else 
		return true;
}
