#include "basetracking.h"

Basetracking::Basetracking(Random *ran, Bfield *bf, Basegeometry *geo)
: fBfield(bf), fRandomgenerator(ran), fGeometry(geo), fIndex(0),
  fLasttime(0), fStarttime(0), fDelta(0.0)
{

}

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

void Basetracking::getB(double time, double* bvec)
{
	if(time > fTracktimes.back())
	{
		throw "Basetracking::getB called with requested time greater than endtime of track";
	}
	else if(time < fTracktimes[0])
	{
		throw "Basetracking::getB called with requested time smaller than starttime of track";
	}

	// this needs the requested times to be in order
	while(time < fTracktimes[fIndex] && fTracktimes[fIndex] < fTracktimes.back())
	{
		fIndex++;
	}
	fDelta = (time - fTracktimes[fIndex-1]) / (fTracktimes[fIndex] - fTracktimes[fIndex-1]);
	Threevector pos = fTrackpositions[fIndex-1] + fDelta * fTrackvelocities[fIndex-1];
	fBfield->eval(time,pos,bvec);
	fLasttime = time;
}

void Basetracking::reset()
{
	fIndex = 0;
	fLasttime = fStarttime;
}

void Basetracking::stepDone()
{
	// this needs the requested times in getB() to be in order
	Threevector lastpos = fTrackpositions[fIndex-1] + fDelta * fTrackvelocities[fIndex-1];
	Threevector lastvel = fTrackvelocities[fIndex-1];	
	fStarttime = fLasttime;
	fIndex = 0;
	fTracktimes.clear();
	fTrackpositions.clear();
	fTrackvelocities.clear();
	fTracktimes.push_back(fLasttime);
	fTrackpositions.push_back(lastpos);
	fTrackvelocities.push_back(lastvel);
}
