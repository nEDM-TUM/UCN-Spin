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
		virtual void initialize();
		virtual Threevector getPosition(double time);
		virtual void makeTrack(double t_start, double h);
		virtual void reset();
		virtual void stepDone(double time);
		
		bool reachedendoftube;
		bool savetrack;
		int Nwallcollisions;
		std::fstream trackparticle;
		std::vector<double> times;
		std::vector<Threevector> positions;
		std::vector<Threevector> axes;
	
	private:
		double vdrift, diffusionconstant, scatteringtime, tend;
		int Nstart;
		bool wasinlastsegment;
		Random *rand;
		Tubegeometry* fTubegeometry;
};
 
#endif
