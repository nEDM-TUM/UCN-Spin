#include "../cylinder.h"
#include "../threevector.h"
#include <iostream>

int main() {
	// how long the simulation should run
	const double T = 300;
	// which stepsize should be used
	const double dt = 1e-3;

	// Cylinder for testing
	Cylinder c(5, 2);

	// Place and velocity for simple test
	Threevector v(2.2, .8, 1);
	Threevector x(1, 1, 1);

	// state for reflect
	bool state[2] = {false, false};

	for (int n = 0; n <= (int) (T/dt); n++) {
		c.reflect(v, x, state);
		
		for (int i = 0; i <= 2; i++)
			x[i] += v[i]*dt;

		std::cout << n*dt << " " << x[0] << " " << x[1] << " " << x[2] << " " << v[0] << " " << v[1] << " " << v[2] << std::endl;
	}
}
