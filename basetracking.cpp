#include "basetracking.h"

/**
 * @class Basetracking
 * Abstract base class, derive from it to create an actual user tracking class.
 * Provides methods to initialize the starting point of a track portion,
 * to reset the tracking, if a timestep was not accepted, 
 * to mark the track as finished and prepare for a new one
 * and to provide the B field along the track
 */

/**
 *
 * @param ran Pointer to the random generator
 * @param bf Pointer to the magnetic field object
 * @param geo Pointer to the geometry object
 */

Basetracking::Basetracking(Random *ran, Basegeometry *geo)
: fRandomgenerator(ran), fGeometry(geo), fIndex(0),
  fLasttime(0), fStarttime(0), fDelta(0.0)
{

}

/**
 * Initializes the track portion. Generates a starting point and sets the
 * starting time to zero
 */

void Basetracking::initialize()
{
	fLasttime = fStarttime = 0.0;
	Threevector pos;
	Threevector vel;
	fGeometry->initialize(vel, pos);
	fTracktimes.push_back(0.0);
	fTrackpositions.push_back(pos);
	fTrackpositions.push_back(vel);
}

/**
 * Returns the magnetic field vector at a certain time
 * @param[in] time The time at which the magnetic field is requested
 * @param[out] bvec Array of 3 doubles that contains the magnetic field vector at the requested time
 */
// TODO: delete
//void Basetracking::getB(double time, double* bvec)
//{
//	if(time > fTracktimes.back())
//	{
//		throw "Basetracking::getB called with requested time greater than endtime of track";
//	}
//	else if(time < fTracktimes[0])
//	{
//		throw "Basetracking::getB called with requested time smaller than starttime of track";
//	}
//
//	// this needs the requested times to be in order
//	while(time < fTracktimes[fIndex] && fTracktimes[fIndex] < fTracktimes.back())
//	{
//		fIndex++;
//	}
//	fDelta = (time - fTracktimes[fIndex-1]) / (fTracktimes[fIndex] - fTracktimes[fIndex-1]);
//	Threevector pos = fTrackpositions[fIndex-1] + fDelta * fTrackvelocities[fIndex-1];
//	fBfield->eval(time,pos,bvec);
//	fLasttime = time;
//}

/**
 * Resets the already constructed track portion to its starting point.
 * This method should be called if the timestep was not accepted.
 */

void Basetracking::reset()
{
	fIndex = 0;
	fLasttime = fStarttime;
}

/**
 * Sets the endpoint and time of the last track portion as starting point and
 * time of the next track portion and deletes the previous track portion.
 * This method should be called if the timestep was accepted.
 */
// TODO: delete
//void Basetracking::stepDone()
//{
//	// this needs the requested times in getB() to be in order
//	Threevector lastpos = fTrackpositions[fIndex-1] + fDelta * fTrackvelocities[fIndex-1];
//	Threevector lastvel = fTrackvelocities[fIndex-1];	
//	fStarttime = fLasttime;
//	fIndex = 0;
//	fTracktimes.clear();
//	fTrackpositions.clear();
//	fTrackvelocities.clear();
//	fTracktimes.push_back(fLasttime);
//	fTrackpositions.push_back(lastpos);
//	fTrackvelocities.push_back(lastvel);
//}
