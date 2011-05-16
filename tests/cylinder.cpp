#include "../cylinder.h"
#include "../threevector.h"
#include "../random.h"
#include <iostream>
#include <cctype>

int main(int argc, char *argv[]) {
	// how long the simulation should run
	const double T = 300;
	// which stepsize should be used
	const double dt = 1e-3;

	// seed for rng
	int seed = 1;
	if (argc == 2)
		seed = atoi(argv[1]);

	Random ran(seed);

	// Cylinder for testing
	Cylinder c(&ran, 5, 2);

	// Place and velocity for simple test
	Threevector v(0, 0, 0);
	Threevector x(0, 0, 0);

	c.initialize(v, x);

	for (int n = 0; n <= (int) (T/dt); n++) {
		if (c.boundsCheck(x))
			c.reflect(v, x);
		
		for (int i = 0; i <= 2; i++)
			x[i] += v[i]*dt;

		std::cout << n*dt << " " << x[0] << " " << x[1] << " " << x[2] << " " << v[0] << " " << v[1] << " " << v[2] << std::endl;
	}
}
