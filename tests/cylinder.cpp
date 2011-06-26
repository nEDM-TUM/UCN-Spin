#include "../cylinder.h"
#include "../fastcylindertracker.h"
#include "../threevector.h"
#include "../random.h"
#include "../debug.h"
#include "../parameters.h"
#include <iostream>
#include <cctype>
#include <unistd.h>

using namespace std;

void outputThreevector(const Threevector& v) {
	for (int i = 0; i < 3; i++)
		cout << v[i];
}

int main(int argc, char *argv[]) {
	// how long the simulation should run
	const double T = 30;

	// how long until the simulation is aborted
	const unsigned int TIMEOUT = 300; // seconds

	alarm(TIMEOUT);

	initialize_debug();

	// seed for rng
	int seed = 1;
	if (argc == 2)
		seed = atoi(argv[1]);

	Random ran(seed);

	Parameters params;
	params.add("CylinderRadius", 0.235);
	params.add("CylinderHeight", 0.12);
	params.add("Temperature", 0.001);
	params.add("ParticleMass", 1.674927e-27);
	params.add("VelocityCutoff", 1.0e100);
	params.add("GravitationConstant", 9.81);
	params.add("DiffusionProbability", 0.);

	// Cylinder for testing
	Cylinder c(params, &ran);

	// Gravitation tracker
	FastCylinderTracker t(params, &ran, &c);
	t.initialize();

	// Make all of the track at once
	t.makeTrack(0, T);

	alarm(TIMEOUT); // Restart timer for writing of output

	// output track
	for (int i = 0; i < t.fTracktimes.size(); i++) {
		cout << t.fTracktimes[i] << " ";
		outputThreevector(t.fTrackpositions[i]);
		cout << " ";
		outputThreevector(t.fTrackvelocities[i]);
		cout << endl;
	}

	return 0;
}
