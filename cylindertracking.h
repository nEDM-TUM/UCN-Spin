#ifndef BASETRACKING_H
#define BASETRACKING_H

#include <vector>
#include "random.h"
#include "bfield.h"
#include "parameters.h"
#include "basegeometry.h"

class CylinderTracking
{
	public:
		CylinderTracking(Random* ran, Bfield* bf, Basegeometry *geo);
		
		virtual void initialize();
		Threevector getPosition(double time);

		/**
		 * Pure virtual function
		 * Creates a track with length @p hmax in time
		 * This method must be provided by a user tracking class
		 */

		virtual void makeTrack(double hmax);
		void reset();
		void stepDone();

	protected:
		Bfield *fBfield;
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
