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
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

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
void setTitles(TH1 &h, const char *x_title, const char *y_title);

int main(int nargs, char** argv)
{
	initialize_debug();

	double firsthtry = 1e-6;
	
	Parameters theParameters;
	theParameters.expectDouble("GyromagneticRatio");
	theParameters.expectInt("NumberOfParticles");
	theParameters.expectDouble("Lifetime");
	theParameters.expectDouble("ErrorGoal");
	theParameters.expectDouble("Flipangle");
	theParameters.expectInt("Seed");
	theParameters.expectDouble("CylinderRadius");
	theParameters.expectDouble("CylinderHeight");
	theParameters.expectDouble("SaveTimeDiff");
	//theParameters.expectDouble("TrackerStepSize");
	theParameters.expectDouble("GravitationConstant");
	//theParameters.expectDouble("CollisionAccuracy");
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

	double approx_b0 = 0; // approximate B0 for expected polarization at end of programm

	// Open output file
	TFile out(generateFileName().c_str(), "new");

	// Set drawing style
	gROOT->SetStyle("Plain");
	gStyle->SetLabelSize(0.025, "xyz");

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
		double P[3] = {-sin(flipangle/180*M_PI), 0.0, cos(flipangle/180*M_PI)};
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
				P[0] = -sin(flipangle/180*M_PI);
				P[1] = 0.0;
				P[2] = cos(flipangle/180*M_PI);
				tracker->initialize();
				approx_b0 += bfield->eval(0)[2];
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

	cout << "Creating histograms..." << endl;

	// Make some histograms
	TH1D end_polarization_hist("end_polarization_hist", "Endpolarisation", 100, 0, 2*M_PI);
	end_polarization.Draw("atan2(polarization.y,polarization.x)+pi>>end_polarization_hist", "", "goff,norm");
	setTitles(end_polarization_hist, "Polarisation [rad]", "");

	TH1D start_velocity_hist("start_velocity_hist", "Anfangsgeschwindigkeit", 100, 0, 1);
	start_velocity_hist.SetBit(TH1::kCanRebin);
	start_tree.Draw("sqrt(velocity.x^2+velocity.y^2+velocity.z^2)>>start_velocity_hist", "", "goff,norm");
	setTitles(start_velocity_hist, "Anfangsgeschwindigkeit [m/s]", "");

	out.Write();

	approx_b0 /= N_particles;
	cout << endl << "Run complete" << endl;
	cout << fixed;
	cout.precision(15);
	cout << "============" << endl;
	cout << "End polarization: (" << end_polarization_hist.GetMean() << " +- " << end_polarization_hist.GetRMS() << ") rad" << endl;
	cout << "Expected:         " << fmod(fmod(-theParameters.getDoubleParam("GyromagneticRatio")*lifetime*approx_b0, 2*M_PI) + 2*M_PI, 2*M_PI) << " rad (approximate value, based on B0 = " << approx_b0 << ")" << endl;
		
	return 0;
}

void setTitles(TH1 &h, const char *x_title, const char *y_title) {
	h.GetXaxis()->SetTitle(x_title);
	h.GetXaxis()->CenterTitle();
	h.GetYaxis()->SetTitle(y_title);
	h.GetYaxis()->CenterTitle();
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
