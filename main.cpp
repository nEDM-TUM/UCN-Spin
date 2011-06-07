#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"

#include "globals.h"
#include "bfield.h"
#include "random.h"
#include "tracking.h"
#include "dopr.h"
#include "derivatives.h"
#include "parameters.h"
#include "gravitationtracker.h"
#include "cylinder.h"
#include "debug.h"

using namespace std;

string generateFileName(void);

int main(int nargs, char** argv)
{
	initialize_debug();

	double firsthtry = 1e-6;
	
	Parameters theParameters;
	theParameters.expectDouble("GyromagneticRatio");
	theParameters.expectInt("NumberOfParticles");
	theParameters.expectDouble("Lifetime");
	theParameters.expectDouble("ErrorGoal");
	theParameters.expectDouble("EfieldMag");
	theParameters.expectDouble("B0Gradient");
	theParameters.expectDouble("GradientOffsetX");
	theParameters.expectDouble("GradientOffsetZ");
	theParameters.expectDouble("B1Gradient");
	theParameters.expectDouble("EDM");
	theParameters.expectDouble("Flipangle");
	theParameters.expectInt("Seed");
	theParameters.expectDouble("CylinderRadius");
	theParameters.expectDouble("CylinderHeight");
	theParameters.expectDouble("B0");
	theParameters.expectDouble("B1");
	theParameters.expectDouble("mu");
	theParameters.expectDouble("radiustube");
	theParameters.expectDouble("vdrift");
	theParameters.expectDouble("sigma");
	theParameters.expectDouble("SolenoidField");
	theParameters.expectDouble("SolenoidCurrent");
	theParameters.expectDouble("SolenoidHeight");
	theParameters.expectDouble("SolenoidRadius");
	theParameters.expectDouble("SaveTimeDiff");

	theParameters.readParameters(cin);

	// TODO: reorder this if possible
	theParameters.add("GyroelectricRatio",theParameters.getDoubleParam("EDM")*0.01*elementarycharge/hbar);

	const double savetimediff = theParameters.getDoubleParam("SaveTimeDiff");
	const int N_particles = theParameters.getIntParam("NumberOfParticles");

	// Open output file
	TFile out(generateFileName().c_str(), "new");

	// Trees for ROOT
	double P_end[3]; // Polarization on end of simulation
	double t_end; // time at end of simulation
	TTree end_polarization("end_polarization", "Polarization at end of run");
	end_polarization.Branch("polarization", P_end, "x/D:y:z");
	end_polarization.Branch("time", &t_end, "t/D");

	#pragma omp parallel
	{
		Random *randgen = new Random(theParameters.getIntParam("Seed"));
		double T = 0.0;
		double flipangle = theParameters.getDoubleParam("Flipangle");
		double P[3] = {0.0,sin(-flipangle/180*M_PI),cos(-flipangle/180*M_PI)};	//the polarization-vector
		double dPdt[3] = {0.0};
		double hdid = 0.0;
		int savetime = 0;

		// TODO
		Cylinder *c = new Cylinder(randgen, 5, 2);
		GravitationTracker *tracker = new GravitationTracker(randgen, c, 9.81); // TODO: parameters

		Bfield *bfield = new Bfield(theParameters,tracker);


		Derivatives * derivatives = new Derivatives(theParameters, bfield);
		double errorgoal = theParameters.getDoubleParam("ErrorGoal");
		Dopr *stepper = new Dopr(firsthtry,	//initial stepsize guess 
					 3,		//dimension of ODE sytem
					 P,
					 dPdt,
					 T,
					 derivatives,
					 errorgoal,	//relative error tolerance
					 errorgoal,	//absolute error tolerance
					 true,		//dense output?
					 tracker
					 );
		
		#pragma omp for
			for (int i = 0; i < N_particles; i++)
			{
				cout << "Particle " << i << "/" << N_particles << endl;
				T = 0.0;
				int Nsteps = 0;
				P[0] = 0.0;
				P[1] = sin(-flipangle/180*M_PI);
				P[2] = cos(-flipangle/180*M_PI);
				tracker->initialize();
				stepper->reset(firsthtry, P, dPdt, T);
				savetime = 0;
				int lifetime = theParameters.getDoubleParam("Lifetime");
				debug << "Will run for " << lifetime << "seconds" << endl;
				while(T <= lifetime) // TODO: Abbruchbedingung
				{
					try
					{
						stepper->step();
						Nsteps++;
						debug << "Step " << Nsteps << " done" << endl;
					}
					catch(char const* error)
					{
						cerr << "ERROR: "<< error << "at t = " << T << endl;
						exit(1);
					}
					
					T = stepper->getT();
					debug << "Stepper TIME: " << T << endl;
					hdid = stepper->getHdid();
					debug << "hdid = " << hdid << endl;

//					while(T >= (st = savetime*savetimediff))
//					{
//						Pol[0] = stepper->dense_out(0,st,hdid);
//						Pol[1] = stepper->dense_out(1,st,hdid);
//						Pol[2] = stepper->dense_out(2,st,hdid);
//						tempTP[4*j] = st;
//						tempTP[4*j+1] = Pol[0];
//						tempTP[4*j+2] = Pol[1];
//						tempTP[4*j+3] = Pol[2];
//						savetime++;
//						j++;
//					}			
				}
				#pragma omp critical
				{
					// fill TTree
					for (int j = 0; j < 3; j++) {
						P_end[j] = stepper->dense_out(j,lifetime,hdid);
					}
					t_end = lifetime;
					end_polarization.Fill();
				}

				cout << "Particle " << (i+1) << ": " << Nsteps << " steps successful, " << stepper->getStepsnottaken() << " steps not taken!" << endl;
			}
		
		delete randgen; randgen = 0;
		delete stepper; stepper = 0;
		delete derivatives; derivatives = 0;
		delete tracker; tracker = 0;
		delete bfield; bfield = 0;
		delete c; c = 0;

	}

	out.Write();
		
	return 0;
}

string generateFileName(void) {
	const size_t MAX_DATE_LEN = 100;
	time_t rawtime;
	char date[MAX_DATE_LEN] = "";

	// get date
	time(&rawtime);
	strftime(date, MAX_DATE_LEN, "%Y-%m-%d-%H-%M-%S-%Z-pid-", gmtime(&rawtime));

	// add pid
	ostringstream o;
	o << "./output/" << date << getpid() << ".root";
	return o.str();
}
