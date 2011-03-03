#include "tracking2.h"
#include "random.h"
#include "bfield.h"
#include <math.h>
#include <iostream>

using namespace std;

int main()
{
	Random *rand = new Random(2545);
	Bfield* B = new Bfield(1e-6, 3e-8, 0.7e-3, 2e-3);
	Tracking *tracker = new Tracking(rand, B, 8.3e-8, 300e-6);
	
	int N = 10000;
	double delta = 2.*1000./(double)N;
	double t = 0.;
	double xyz[3];
	double vxyz[3];
	tracker->initialize();
	for(int i=1; i<=N; i++)
	{
		tracker->getPosAndVel(t, xyz, vxyz);
		cout << sqrt(vxyz[0]*vxyz[0]+vxyz[1]*vxyz[1]+vxyz[2]*vxyz[2]) << endl;
		t += delta*rand->uniform();
	}
	delete rand;
	delete B;
	delete tracker;
	return 0;
}
