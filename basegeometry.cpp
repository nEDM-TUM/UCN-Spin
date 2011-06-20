#include <basegeometry.h>
#include <cmath>

void diffuse(Threevector &v, const Threevector &x) {
	double vx, vy, vz; // New velocity
	
	#pragma omp critical
	{
		vx = fRandom->gaussian(1);
		vy = fRandom->gaussian(1);
		vz = fRandom->gaussian(1);
	}

	// factor which has to be multiplied to vx, vy, vz to conserve velocity
	const double v_fact = sqrt(v.magsquare() / (vx*vx + vy*vy + vz*vz));

	return Threevector(v_fact*vx, v_fact*vy, v_fact*vz);
}
