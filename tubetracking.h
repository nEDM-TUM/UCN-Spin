#ifndef TUBETRACKING_H
#define TUBETRACKING_H
#include <vector>
#include <fstream>
#include <iostream>
#include "random.h"
#include "parameters.h"
#include "basetracking.h"
#include "threevector.h"
#include "tubegeometry.h"


class Tubetracking : public Basetracking {	
	public:
		Tubetracking(Random*, Tubegeometry * geo, Parameters& theParameters);
		~Tubetracking();
		
// Ruft Initialize von Tubegeometry auf. Fügt den 1. axialen Vektor zu 'axes' 
// und den Startpunkt zu 'positions' hinzu. 
		virtual void initialize();

// Gibt die Position des Xe-Atoms zu einer bestimmten Zeit 'time' zurück.
		virtual Threevector getPosition(double time);
		virtual void makeTrack(double t_start, double h);
		virtual void reset();

// Löscht den Inhalt aus den Vektoren 'positions', 'axes', 'times', 'roots' und fügt 
// gleichzeitig die jeweils letzten Werte wieder hinzu.
		virtual void stepDone(double time);
		
		bool reachedendoftube;
		bool savetrack;
		std::fstream trackparticle;
// In den Vektoren sollen die bisherigen ausgewürfelten Werte, die zu einer "guten"
// Position geführt haben, gespeichert werden. 
		std::vector<double> times;
		//std::vector<double> roots;
		std::vector<Threevector> positions;
		std::vector<Threevector> axes;
	
	private:
		double vdrift, mu, sigma, scatteringtime, tend;
		int Nstart;
		bool wasinlastsegment;
		Random *rand;
		Tubegeometry* fTubegeometry;
};
 
#endif
