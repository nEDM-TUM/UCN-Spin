using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include "globals.h"
#include "bfield.h"
#include "random.h"
#include "tracking.h"
#include "dopr.h"
#include "derivatives.h"
#include "parameters.h"


int main(int nargs, char** argv)
{
	double savetimediff = 2e-2, firsthtry = 1e-6;
	int iP=0;
	double *TP_shared = NULL, *P_end = NULL;
	
	Parameters theParameters;
	std::vector<std::string> parameternames;
	parameternames.push_back("NumberOfParticles");
	parameternames.push_back("Lifetime");
	parameternames.push_back("ErrorGoal");
	parameternames.push_back("EfieldMag");
	parameternames.push_back("B0Gradient");
	parameternames.push_back("GradientOffsetX");
	parameternames.push_back("GradientOffsetZ");
	parameternames.push_back("B1Gradient");
	parameternames.push_back("EDM");
	parameternames.push_back("Flipangle");
	parameternames.push_back("Seed");

	if(nargs == parameternames.size()+1)
	{
		for(size_t i=1; i<= parameternames.size(); i++)
		{
			if(parameternames.at(i-1)=="NumberOfParticles" || parameternames.at(i-1)=="Seed")	//for integer parameters
			{
				theParameters.add(parameternames.at(i-1),atoi(argv[i]));
			}
			else
			{
				theParameters.add(parameternames.at(i-1),atof(argv[i]));
			}
		}
	}
	else
	{
		cout << "wrong number of arguments passed: " << nargs << endl;
		exit(1);
	}

	theParameters.add("GyroelectricRatio",theParameters.getDoubleParam("EDM")*0.01*elementarycharge/hbar);
	theParameters.add("Radius",3.0e-4);
	theParameters.add("B0",1.0e-6);
	theParameters.add("B1",1.0e-9);

	const int J = 4*((int)(theParameters.getDoubleParam("Lifetime")/savetimediff) + 14);
	int N_particles = theParameters.getIntParam("NumberOfParticles");
	TP_shared = new double[J];
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
		
	Random *randgen = new Random(theParameters.getIntParam("Seed"));

	#pragma omp parallel
	{
		double T = 0.0;
		double Told = 0.0;
		double flipangle = theParameters.getDoubleParam("Flipangle");
		double P[3] = {0.0,sin(-flipangle/180*M_PI),cos(-flipangle/180*M_PI)};	//the polarization-vector
		double dPdt[3] = {0.0};
		double hdid = 0.0;
		int savetime = 0;
		double *tempTP = new double[J];
		double *TP = new double[J];
		for(int z=0; z<J; z++)
		{
			tempTP[z] = 0.0;
			TP[z] = 0.0;
		}
		int i=0,j=0;
		double Pol[3] = {0.0,0.0,0.0};

		for(int z=0; z<J; z++)
		{
			TP[z] = 0.0;
		}
		
		string rotationBfilename = "./RotatingBfield.txt";	//File which contains B-field frequencies and times
		Bfield *bfield = new Bfield(rotationBfilename,theParameters);
													
		Tracking *tracker = new Tracking(randgen, 		//Pointer to the random-generator
						 bfield,		//Pointer to the magnetic field
						 theParameters);

		Derivatives * derivatives = new Derivatives(tracker);
		double errorgoal = theParameters.getDoubleParam("ErrorGoal");
		Dopr *stepper = new Dopr(firsthtry,	//initial stepsize guess 
					 3,		//dimension of ODE sytem
					 P,
					 dPdt,
					 T,
					 derivatives,
					 errorgoal,	//relative error tolerance
					 errorgoal,	//absolute error tolerance
					 true);		//dense output?
		
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
				stepper->reset(firsthtry, P, dPdt, T);
				savetime = 0;
				j = 0;
				int lifetime1 = theParameters.getDoubleParam("Lifetime");
				while(T < lifetime1)
				{
					try
					{
						stepper->step();
						Nsteps++;
					}
					catch(char const* error)
					{
						cerr << "ERROR: "<< error << "at t = " << T << endl;
						exit(1);
					}
					Told = T;
					T = stepper->getT();			
					hdid = stepper->getHdid();
					double st = 0.0;
					while(T >= (st = savetime*savetimediff))
					{
						Pol[0] = stepper->dense_out(0,st,hdid);
						Pol[1] = stepper->dense_out(1,st,hdid);
						Pol[2] = stepper->dense_out(2,st,hdid);
						tempTP[4*j] = st;
						tempTP[4*j+1] = Pol[0];
						tempTP[4*j+2] = Pol[1];
						tempTP[4*j+3] = Pol[2];
						savetime++;
						j++;
					}			
				}
				#pragma omp critical
				{
					P_end[3*iP]   = stepper->dense_out(0,lifetime1,hdid);
					P_end[3*iP+1] = stepper->dense_out(1,lifetime1,hdid);
					P_end[3*iP+2] = stepper->dense_out(2,lifetime1,hdid);
					iP++;
				}
				for(int z=0; z<J; z++)
				{
					TP[z] += tempTP[z];
				}
				cout << Nsteps << " steps successful, " << stepper->getStepsnottaken() << " steps not taken!" << endl;
			}
		#pragma omp critical
		{
			for(int z=0; z<J; z++)
			{
				TP_shared[z] += TP[z];
			}
		}
		
		delete stepper; stepper = NULL;
		delete derivatives; derivatives = NULL;
		delete tracker; tracker = NULL;
		delete bfield; bfield = NULL;
		
		delete[] tempTP; tempTP = NULL;
		delete[] TP; TP = NULL;
	}
	string ofilename = "";
	string ofilename1 = "";
	stringstream NumberString;
	NumberString << N_particles;
	ofilename1 = "("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("Lifetime");
	ofilename = "Lifetime("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("ErrorGoal");
	ofilename += "Err("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("EfieldMag");
	ofilename += "Efield("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("B0Gradient");
	ofilename += "Gradient("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("GradientOffsetX");
	ofilename += "x-offset("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("GradientOffsetZ");
	ofilename += "z-offset("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("B1Gradient");
	ofilename += "B1-gradient("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("B1");
	ofilename += "B1("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("EDM");
	ofilename += "EDM("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getDoubleParam("Flipangle");
	ofilename += "Flip("+NumberString.str()+")_";
	NumberString.str("");
	NumberString.clear();
	NumberString << theParameters.getIntParam("Seed");
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
	
	ofstream Pend(("./output/Pend/Pend"+ofilename1+ofilename).c_str());
	Pend.precision(15);
	for(int z=0; z<N_particles; z++)
	{
		Pend << P_end[3*z] << "\t" << P_end[3*z+1] << "\t" << P_end[3*z+2] << endl;
	}
	Pend.close();
	
	delete[] TP_shared; TP_shared = NULL;
	delete[] P_end; P_end = NULL;
	return 0;
}
