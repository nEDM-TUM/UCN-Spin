using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include "globals.h"
#include "bfield.h"
#include "random.h"
#include "tracking.h"
#include "spintracking.h"
#include "derivatives.h"
#include "dipolefield.h"


int main(int nargs, char** argv)
{
	const double hbar = 1.054571682364455e-34;
	const double elementarycharge = 1.60218e-19;
	double savetimediff = 2e-2, Radius = 3e-4, error = 0.0, firsthtry = 1e-6, lifetime1 = 0.0, cylindercurrent = 0.0;
	double g = 0.0, B1_g = 0.0, E1 = 0., B1 = 1e-9, B0 = 1e-6, D = 0.0, EDM = 0.0, flipangle = 0.0, gyroelect = 0.0;
	double xyz0[3] = {0.0,0.0,0.0};
	double Distances[100] = {0.0};
	bool rotatedipole;
	unsigned long int seed = 123;
	int iP=0;
	int Ndistances = 0;
	int N_particles = 2;
	double *TP_shared = NULL, **BDip_shared = NULL, *P_end = NULL;
	ifstream Distancefile;
	Distancefile.open("./Distances.txt");	// file contains the distances in x-direction, where the signal is measured
	if(!Distancefile)
	{
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}
	while(!(Distancefile.eof()))
	{
		Distancefile >> Distances[Ndistances];
		Ndistances++;
	}
	Ndistances--;
	
	if(nargs == 15)								// COMMANDLINE ARGUMENTS:
	{
		N_particles = atoi(argv[1]);			// Number of particles
		lifetime1 = atof(argv[2]);				// Lifetime for each particle
		error = atof(argv[3]);					// RK error
		cylindercurrent = atof(argv[4]);		// current in the cylinder induced by rotating efield
		E1 = atof(argv[5]);						// magnitude of efield (V/m)
		g = atof(argv[6]);						// gradient in z-direction (T/m)
		xyz0[0] = atof(argv[7]);				// gradient offset in x-direction (m)
		xyz0[2] = atof(argv[8]);				// gradient offset in z-direction (m)
		B1_g = atof(argv[9]);					// B1-gradient (T/m)
		D = atof(argv[10]);						// self diffusion coefficient
		EDM = atof(argv[11]);					// Electric dipole moment (e cm)
		flipangle = atof(argv[12]);				// The flipangle (degree)
		rotatedipole = atoi(argv[13]);			// Rotate the dipolefield through an EDM? 
												   // (if 1, the simulated droplet will not be exposed to
												   // an EDM-rotation, only the neighbour-droplet will be 
												   // rotated. if 0 then its the other way)
		seed = atol(argv[14]);					// RNG seed
	}
	else
	{
		cout << "wrong number of arguments passed: " << nargs << endl;
		exit(1);
	}
	const int J = 4*((int)(lifetime1/savetimediff) + 14);
	gyroelect = EDM*0.01*elementarycharge/hbar;		// now in SI units (m/(V s))
	TP_shared = new double[J];
	BDip_shared = new double*[Ndistances];
	BDip_shared[0] = new double[J*Ndistances];
	for(int z=1; z<Ndistances; z++)
	{
		BDip_shared[z] = BDip_shared[z-1] + J;
	}
	P_end = new double[3*N_particles];
	#pragma omp parallel for
	for(int z=0; z<J; z++)
	{
		TP_shared[z] = 0.0;
	}
	#pragma omp parallel for
	for(int z=0; z<3*N_particles; z++)
	{
		P_end[z] = 0.0;
	}
	for(int y=0; y<Ndistances; y++)
	{
		#pragma omp parallel for
		for(int z=0; z<J; z++)
		{
			BDip_shared[y][z] = 0.0;
		}
	}	
		
	Random *randgen = new Random(seed);

	#pragma omp parallel
	{
		double T = 0.0;
		double Told = 0.0;
		double P[3] = {0.0,sin(-flipangle/180*M_PI),cos(-flipangle/180*M_PI)};	//the polarization-vector
		double hdid = 0.0;
		int savetime = 0;
		double *tempTP = new double[J];
		double *TP = new double[J];
		for(int z=0; z<J; z++)
		{
			tempTP[z] = 0.0;
			TP[z] = 0.0;
		}
		double **tempBDip = new double*[Ndistances];
		double **BDip = new double*[Ndistances];
		tempBDip[0] = new double[J*Ndistances];
		BDip[0] = new double[J*Ndistances];
		for(int z=1; z<Ndistances; z++)
		{
			tempBDip[z] = tempBDip[z-1] + J;
			BDip[z] = BDip[z-1] + J;
		}
		for(int y=0; y<Ndistances; y++)
		{
			for(int z=0; z<J; z++)
			{
				tempBDip[y][z] = 0.0;
				BDip[y][z] = 0.0;
			}
		}
		int i=0,j=0;
		double Pol[3] = {0.0,0.0,0.0};
		double Bf[3] = {0.0,0.0,0.0};

		for(int z=0; z<J; z++)
		{
			TP[z] = 0.0;
		}
		
		string rotationBfilename = "./RotatingBfield.txt";	//File which contains B-field frequencies and times
		string rotationDipfilename = "./RotatingDipolefield.txt";
		string dipoleLocations = "./DipoleLocations.txt";
		Bfield *bfield = new Bfield(rotationBfilename,
									rotationDipfilename,
									dipoleLocations,
									B0,					//B0-holding-field in (T)
									B1,					//Magnitude of B1-field in (T)
									cylindercurrent,	//Cylindrical current in (A)
									0.75e-3,			//Radius of the current-cylinder in (m)
									0.002,				//Full height of the current-cylinder in (m)
									E1,					//Strength of the E-field (V/m)
									g,					//B0-gradient in z-direction (T/m)
									xyz0,				//B0-gradient offset (m,m,m)
									B1_g,				//B1-gradient (T/m)
									gyroelect,			//GYROELECTRIC RATIO
									flipangle/180*M_PI,	//the flipangle of the (neighbour-droplet) dipolefield
									rotatedipole);		//wheter an EDM is applied to the simulated or the neighbour-droplet
													
		Tracking *tracker = new Tracking(randgen, 		//Pointer to the random-generator
										 bfield,		//Pointer to the magnetic field
										 D,				//Self-Diffusion-constant of Xe (8.3e-8)
										 Radius);		//Radius of the Xe-droplet

		Derivatives * derivatives = new Derivatives(tracker);
		Spintracking *spintracker = new Spintracking(P, T, error, firsthtry, derivatives, true);
		Dipolefield *Dipole = new Dipolefield(B0);
		
		#pragma omp for
			for(i=0; i<N_particles; i++)
			{
				cout << "Particle " << i << "/" << N_particles << endl;
				T = 0.0;
				Told = 0.0;
				int Nsteps = 0;
				P[0] = 0.0;
				P[1] = sin(-flipangle/180*M_PI);
				P[2] = cos(-flipangle/180*M_PI);
				tracker->initialize();
				spintracker->reset(P, T, firsthtry);
				savetime = 0;
				j = 0;
				while(T < lifetime1)
				{
					try
					{
						spintracker->precess();
						Nsteps++;
					}
					catch(char const* error)
					{
						cerr << "ERROR: "<< error << "at t = " << T << endl;
						exit(1);
					}
					Told = T;
					T = spintracker->getStepperTime();			
					hdid = spintracker->getHdid();
					double st = 0.0;
					while(T >= (st = savetime*savetimediff))
					{
						Pol[0] = spintracker->dense_out(0,st,hdid);
						Pol[1] = spintracker->dense_out(1,st,hdid);
						Pol[2] = spintracker->dense_out(2,st,hdid);
						tempTP[4*j] = st;
						tempTP[4*j+1] = Pol[0];
						tempTP[4*j+2] = Pol[1];
						tempTP[4*j+3] = Pol[2];
						double *r0 = tracker->getLastPos();
						for(int z=0; z<Ndistances; z++)
						{
							double r[3]	= {Distances[z]-r0[0],r0[1],r0[2]};
							Dipole->getField(st,Bf,r,Pol);
							tempBDip[z][4*j] = st;
							tempBDip[z][4*j+1] = Bf[0];
							tempBDip[z][4*j+2] = Bf[1];
							tempBDip[z][4*j+3] = Bf[2];
						}
						savetime++;
						j++;
					}			
				}
				#pragma omp critical
				{
					P_end[3*iP]   = spintracker->dense_out(0,lifetime1,hdid);
					P_end[3*iP+1] = spintracker->dense_out(1,lifetime1,hdid);
					P_end[3*iP+2] = spintracker->dense_out(2,lifetime1,hdid);
					iP++;
				}
				for(int z=0; z<J; z++)
				{
					TP[z] += tempTP[z];
					for(int y=0; y<Ndistances; y++)
					{
						BDip[y][z] += tempBDip[y][z];
					}
				}
				cout << Nsteps << " steps successful, " << spintracker->getStepsnottaken() << " steps not taken!" << endl;
			}
		#pragma omp critical
		{
			for(int z=0; z<J; z++)
			{
				TP_shared[z] += TP[z];
				for(int y=0; y<Ndistances; y++)
				{
					BDip_shared[y][z] += BDip[y][z];
				}
			}
		}
		
		delete spintracker; spintracker = NULL;
		delete derivatives; derivatives = NULL;
		delete tracker; tracker = NULL;
		delete bfield; bfield = NULL;
		delete Dipole; Dipole = NULL;
		
		delete[] tempTP; tempTP = NULL;
		delete[] TP; TP = NULL;
		delete[] tempBDip[0]; tempBDip[0] = NULL;
		delete[] tempBDip; tempBDip = NULL;
		delete[] BDip; BDip = NULL;
	}
	string ofilename = "";
	string ofilename1 = "";
	stringstream NumberString;
	NumberString << N_particles;
	ofilename1 = "("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << lifetime1;
	ofilename = "Lifetime("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << error;
	ofilename += "Err("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << cylindercurrent;
	ofilename += "Current("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << E1;
	ofilename += "Efield("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << g;
	ofilename += "Gradient("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << xyz0[0];
	ofilename += "x-offset("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << xyz0[2];
	ofilename += "z-offset("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << B1_g;
	ofilename += "B1-gradient("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << B1;
	ofilename += "B1("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << D;
	ofilename += "Diffusion("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << EDM;
	ofilename += "EDM("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << flipangle;
	ofilename += "Flip("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << rotatedipole;
	ofilename += "Rotatedip("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << seed;
	ofilename += "seed("+NumberString.str()+").dat";
	NumberString.str("");
	NumberString.clear();

	ofstream Poutput(("./output/Pol/Pol"+ofilename1+ofilename).c_str());
	Poutput.precision(15);
	for(int z=0; z<(J-1)/4; z++)
	{
		Poutput << TP_shared[z*4]/N_particles << "\t" << TP_shared[z*4+1]/N_particles << "\t" << TP_shared[z*4+2]/N_particles 
				<< "\t" << TP_shared[z*4+3]/N_particles << endl;
	}
	Poutput.close();
	
	for(int y=0; y<Ndistances; y++)
	{
		NumberString << Distances[y];
		ofstream BDipoutput(("./output/Dip/Dip"+ofilename1+"Distance("+NumberString.str()+")_"+ofilename).c_str());
		BDipoutput.precision(15);
		for(int z=0; z<(J-1)/4; z++)
		{
			BDipoutput << BDip_shared[y][z*4]/N_particles << "\t" << BDip_shared[y][z*4+1]/N_particles << "\t" 
			<< BDip_shared[y][z*4+2]/N_particles << "\t" << BDip_shared[y][z*4+3]/N_particles << endl;
		}
		BDipoutput.close();
		NumberString.str("");
		NumberString.clear();
	}
	
	ofstream Pend(("./output/Pend/Pend"+ofilename1+ofilename).c_str());
	Pend.precision(15);
	for(int z=0; z<N_particles; z++)
	{
		Pend << P_end[3*z] << "\t" << P_end[3*z+1] << "\t" << P_end[3*z+2] << endl;
	}
	Pend.close();
	
	delete[] TP_shared; TP_shared = NULL;
	delete[] BDip_shared[0]; BDip_shared[0] = NULL;
	delete[] BDip_shared; BDip_shared = NULL;
	delete[] P_end; P_end = NULL;
	return 0;
}
