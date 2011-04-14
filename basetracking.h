#ifndef BASETRACKING_H
#define BASETRACKING_H

#include <vector>
#include "random.h"
#include "parameters.h"
#include "basegeometry.h"
#include "threevector.h"

class Basetracking
{
	public:
		Basetracking(Random* ran, Basegeometry *geo);
		
		virtual void initialize();
		virtual Threevector getPosition(double time) = 0;

		/**
		 * Pure virtual function
		 * Creates a track with length @p hmax in time
		 * This method must be provided by a user tracking class
		 */

		virtual void makeTrack(double hmax) = 0;
		void reset();
		void stepDone();

		std::vector<double> fTracktimes;
		std::vector<Threevector> fTrackpositions;
		std::vector<Threevector> fTrackvelocities;

	protected:
		Random *fRandomgenerator;
		Basegeometry *fGeometry;
		int fIndex;
		double fLasttime;
		double fStarttime;
};
#endif 
