#include "basetracking.h"

/**
 * @class Basetracking
 * Abstract base class, derive from it to create an actual user tracking class.
 * Provides methods to initialize the starting point of a track portion,
 * to reset the tracking, if a timestep was not accepted and
 * to mark the track as finished and prepare for a new one.
 */

/**
 *
 * @param ran Pointer to the random generator
 * @param geo Pointer to the geometry object
 */

Basetracking::Basetracking(Random *ran, Basegeometry *geo)
: fRandomgenerator(ran), fGeometry(geo), fIndex(0),
  fLasttime(0), fStarttime(0)
{
}

Basetracking::~Basetracking()
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
	// TODO: how should vel and pos be initialized???
	fGeometry->initialize(vel, pos);
	fTracktimes.push_back(0.0);
	fTrackpositions.push_back(pos);
	fTrackpositions.push_back(vel);
}

/**
 * Resets the already constructed track portion to its starting point.
 * This method should be called if the timestep was not accepted.
 */
void Basetracking::reset()
{
	fIndex = 0;
	fLasttime = fStarttime;
}
