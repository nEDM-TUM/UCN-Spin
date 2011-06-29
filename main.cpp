#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <limits>
#include <unistd.h>
#include <cassert>

#include "TFile.h"
#include "TTree.h"

#include "globals.h"
#include "superpositionfield.h"
#include "random.h"
#include "tracking.h"
#include "dopr.h"
#include "derivatives.h"
#include "parameters.h"
//#include "gravitationtracker.h"
#include "fastcylindertracker.h"
#include "cylinder.h"
#include "debug.h"
#include "exceptions.h"
#include "timeout.h"

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
	theParameters.expectDouble("TrackerStepSize");
	theParameters.expectDouble("GravitationConstant");
	theParameters.expectDouble("CollisionAccuracy");
	theParameters.expectInt("Timeout");
	theParameters.expectInt("ParticleDataNum");
	theParameters.expectDouble("Temperature");
	theParameters.expectDouble("ParticleMass");
	theParameters.expectDouble("VelocityCutoff");
	theParameters.expectDouble("DiffusionProbability");
	theParameters.expectDouble("MinDiffusionAngle");


	theParameters.readParameters(cin);

	const int N_particles = theParameters.getIntParam("NumberOfParticles");
	const int save_every = N_particles / theParameters.getIntParam("ParticleDataNum") + 1;
	const double savetimediff = theParameters.getDoubleParam("SaveTimeDiff");
	const double lifetime = theParameters.getDoubleParam("Lifetime");

	// Open output file
	TFile out(generateFileName().c_str(), "new");

	// Trees for ROOT
	double P_end[3]; // Polarization on end of simulation
	double t_end; // time at end of simulation
	TTree end_polarization("end_polarization", "Polarization at end of run");
	end_polarization.Branch("polarization", P_end, "x/D:y:z");
	end_polarization.Branch("time", &t_end, "t/D");

	double random_val;
	TTree random_tree("random", "Random numbers from start of each run");
	random_tree.Branch("random", &random_val, "random/D");

	float start_velocity[3];
	float start_position[3];
	TTree start_tree("start_values", "Place and velocity at start of run");
	start_tree.Branch("velocity", start_velocity, "x:y:z");
	start_tree.Branch("position", start_position, "x:y:z");

	UInt_t particle_number; // 32 bit unsigned integer
	float particle_time;
	float particle_position[3];
	float particle_velocity[3];
	float particle_polarization[3];
	float particle_field[3];
	TTree particle_tree("particle_data", "Time resolved data for particles");
	particle_tree.Branch("number", &particle_number, "n/i");
	particle_tree.Branch("time", &particle_time, "t");
	particle_tree.Branch("position", particle_position, "x:y:z");
	particle_tree.Branch("velocity", particle_velocity, "x:y:z");
	particle_tree.Branch("polarization", particle_polarization, "x:y:z");
	particle_tree.Branch("bfield", particle_field, "x:y:z");

	Random seed_generator(theParameters.getIntParam("Seed"));

	#pragma omp parallel
	{
		Random *randgen;
		#pragma omp critical
		{
			randgen = new Random(seed_generator.generate_seed());
		}
		double T = 0.0;
		double flipangle = theParameters.getDoubleParam("Flipangle");
		double P[3] = {0.0,sin(-flipangle/180*M_PI),cos(-flipangle/180*M_PI)};	//the polarization-vector
		double dPdt[3] = {0.0};
		double hdid = 0.0;

		Cylinder *c = new Cylinder(theParameters, randgen);
		//GravitationTracker *tracker = new GravitationTracker(theParameters, randgen, c); // TODO
		FastCylinderTracker *tracker = new FastCylinderTracker(theParameters, randgen, c);

		Bfield *bfield = new SuperpositionField(tracker, std::string("fields.dat"));


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
			// save polarization during run because
			// dense_out cannot be used later
			float savetime = 0;
			vector<float> temp_polarization_t;
			vector<float> temp_polarization_x;
			vector<float> temp_polarization_y;
			vector<float> temp_polarization_z;
			if (i % save_every == 0) {
				temp_polarization_t.reserve(1.1*lifetime/savetimediff);
				temp_polarization_x.reserve(1.1*lifetime/savetimediff);
				temp_polarization_y.reserve(1.1*lifetime/savetimediff);
				temp_polarization_z.reserve(1.1*lifetime/savetimediff);
			}

			try {
				Timeout timeout(theParameters.getIntParam("Timeout"));

				#pragma omp critical
				{
					cout << "Particle " << i << "/" << N_particles << endl;
				}
				#pragma omp critical
				{
					random_val = randgen->uniform_double(0, 1);
					random_tree.Fill();
				}
				T = 0.0;
				int Nsteps = 0;
				P[0] = 0.0;
				P[1] = sin(-flipangle/180*M_PI);
				P[2] = cos(-flipangle/180*M_PI);
				tracker->initialize();
				#pragma omp critical
				{
					debug << "Got from initialize: x = " << tracker->fTrackpositions[0].toString() << endl;
					debug << "Got from initialize: v = " << tracker->fTrackvelocities[0].toString() << endl;

					for (int i = 0; i < 3; i++) {
						assert(tracker->fTrackpositions.size() == 1);
						start_position[i] = tracker->fTrackpositions[0][i];
						assert(tracker->fTrackvelocities.size() == 1);
						start_velocity[i] = tracker->fTrackvelocities[0][i];
					}

					start_tree.Fill();
				}
				stepper->reset(firsthtry, P, dPdt, T);
				savetime = 0;
				debug << "Will run for " << lifetime << "seconds" << endl;
				while(T <= lifetime)
				{
					timeout.check();

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

					// save time resolved polarization to temporary storage
					// until end of run.
					while ((i % save_every == 0) && savetime <= T) {
						temp_polarization_t.push_back(savetime);
						temp_polarization_x.push_back(stepper->dense_out(0, savetime, hdid));
						temp_polarization_y.push_back(stepper->dense_out(1, savetime, hdid));
						temp_polarization_z.push_back(stepper->dense_out(2, savetime, hdid));
						savetime += savetimediff;
					}
				} // while

				// write end polarization
				#pragma omp critical
				{
					// fill TTree
					for (int j = 0; j < 3; j++) {
							P_end[j] = stepper->dense_out(j,lifetime,hdid);
					}
					t_end = lifetime;
					end_polarization.Fill();
				}

				// write time resolved data
				if (i % save_every == 0) {
					#pragma omp critical
					{
						assert(temp_polarization_x.size() == temp_polarization_y.size() && temp_polarization_y.size() == temp_polarization_z.size());
						assert(temp_polarization_z.size() == temp_polarization_t.size());

						particle_number = i; // i of particle loop
						for (unsigned int j = 0; j < temp_polarization_t.size(); j++) {
							particle_time = temp_polarization_t[j];
							particle_polarization[0] = temp_polarization_x[j];
							particle_polarization[1] = temp_polarization_y[j];
							particle_polarization[2] = temp_polarization_z[j];
							
							const Threevector pos   = tracker->getPosition(particle_time);
							const Threevector vel   = tracker->getVelocity(particle_time);
							const Threevector field = bfield->eval(particle_time);
							for (int k = 0; k < 3; k++) {
								particle_position[k] = pos[k];
								particle_velocity[k] = vel[k];
								particle_field[k] = field[k];
							}

							particle_tree.Fill();
						}
						temp_polarization_t.clear();
						temp_polarization_x.clear();
						temp_polarization_y.clear();
						temp_polarization_z.clear();
					}
				}

				#pragma omp critical
				{
					cout << "Particle " << (i+1) << ": " << Nsteps << " steps successful, " << stepper->getStepsnottaken() << " steps not taken!" << endl;
				}

			} // try
			catch (const Exception& e) {
				#pragma omp critical
				{
					cout << "Exception for particle " << i << ": " << e.what() << endl;
				}
			}

		} // for

		
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
