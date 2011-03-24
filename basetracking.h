#ifndef BASETRACKING_H
#define BASETRACKING_H

#include <vector>
#include "random.h"
#include "bfield.h"
#include "parameters.h"
#include "basegeometry.h"

class Basetracking
{
	public:
		Basetracking(Random* ran, Basegeometry *geo);
		
		virtual void initialize();
		Threevector getPosition(double time);

		/**
		 * Pure virtual function
		 * Creates a track with length @p hmax in time
		 * This method must be provided by a user tracking class
		 */

		virtual void makeTrack(double hmax) = 0;
		void reset();
		void stepDone();

	protected:
		Random *fRandomgenerator;
		Basegeometry *fGeometry;
		vector<double> fTracktimes;
		vector<Threevector> fTrackpositions;
		vector<Threevector> fTrackvelocities;
		int fIndex;
		double fLasttime;
		double fStarttime;
		double fDelta;
};
#endif 
