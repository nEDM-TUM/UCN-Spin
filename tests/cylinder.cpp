#include "../cylinder.h"
#include "../threevector.h"
#include "../random.h"
#include "../gravitationtracker.h"
#include "../debug.h"
#include <iostream>
#include <cctype>

using namespace std;

inline void outputThreevector(const Threevector &t) {
	cout << t[0] << " " << t[1] << " " << t[2];
}

int main(int argc, char *argv[]) {
	// how long the simulation should run
	const double T = 300;

	initialize_debug();

	// seed for rng
	int seed = 1;
	if (argc == 2)
		seed = atoi(argv[1]);

	Random ran(seed);

	// Cylinder for testing
	Cylinder c(&ran, 5, 2);

	// Gravitation tracker
	GravitationTracker t(&ran, &c, 9.81);
	t.initialize();

	// Make all of the track at once
	t.makeTrack(0, T);

	// output track
	for (int i = 0; i < t.fTracktimes.size(); i++) {
		cout << t.fTracktimes[i] << " ";
		outputThreevector(t.fTrackpositions[i]);
		cout << " ";
		outputThreevector(t.fTrackvelocities[i]);
		cout << endl;
	}
}
