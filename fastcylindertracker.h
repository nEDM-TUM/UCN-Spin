#ifndef _FASTCYLINDERTRACKER_H
#define _FASTCYLINDERTRACKER_H

#include <vector>
#include <utility>
#include "basetracking.h"
#include "cylinder.h"

class FastCylinderTracker : public Basetracking
{
	public:
		FastCylinderTracker(const Parameters &params, Random* ran, Cylinder *geo);
		virtual ~FastCylinderTracker();

		virtual void initialize();
		Threevector getPosition(double time);
		Threevector getVelocity(double time);
		void makeTrack(double t_start, double h);

	private:
		enum Surface {None, Top, Bottom, Radius};
		Cylinder *fGeometry;
		const double H;  //< height of the cylinder
		const double R;  //< radius of the cylinder
		const double R2; //< squared radius of the cylinder
		const double g;  //< acceleration due to gravitation

		Surface fLastCollisionSurface;
		const double fDiffuseProbability;

		void addSolutions(std::vector< std::pair<double, Surface> > &v, const double a, const double b, const double c,
				const double t0, const Surface current) const;

		/**
		 * Find the index in fTracktimes for @time.
		 * Returns i with <tt>fTracktimes[i] <= time <= fTrackpositions[i+1]</tt>.
		 */
		unsigned int findIndex(const double time);
};
#endif // _FASTCYLINDERTRACKER_H
