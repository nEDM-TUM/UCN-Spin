#include "tubetracking.h"
#include "debug.h"
#include <cmath>


Tubetracking::Tubetracking(Random *ran, Tubegeometry * geo, Parameters& theParameters)
:Basetracking(ran, geo), rand(ran), fTubegeometry(geo)
{
	vdrift = theParameters.getDoubleParam("vdrift");
	diffusionconstant = theParameters.getDoubleParam("diffusionconstant");
	scatteringtime = 0;
	Nstart = 0;
	Nwallcollisions = 0;
	wasinlastsegment = false;
	reachedendoftube = false;
	tend = 0;
	trackparticle.open("track.txt");
}

Tubetracking::~Tubetracking(){
	trackparticle.close();
}

void Tubetracking::initialize(){
    positions.clear();
    axes.clear();
    times.clear();
    wasinlastsegment = false;
    reachedendoftube = false;
    tend = 0;
    Nstart = 0;
    Nwallcollisions = 0;
	Threevector v,x;
	v = Threevector();
	x = Threevector(); 
	fTubegeometry->initialize(v,x); 
	x = Threevector(0.,0.,0.);
	positions.push_back(x);
	axes.push_back(v);
	times.push_back(0.0);
}


Threevector Tubetracking::getPosition(double time){
	Threevector pos;	
	
	if (time < 1e-50){
		return positions.back();
	}
	
	while (times.back() < time && wasinlastsegment == false) {
		double tnew, sigma, scatteringlength_1, scatteringlength_2, scatteringlength_3;
		Threevector positionnew, axisnew, control, scatteringvector;
		tnew = times.back() + scatteringtime; 
		control = Threevector (0,0,0);
		sigma = sqrt(2*diffusionconstant*scatteringtime);
		
		//scatteringlength_1 = fRandomgenerator->gaussian(sigma); 
		//scatteringlength_2 = fRandomgenerator->gaussian(sigma);
		//scatteringlength_3 = fRandomgenerator->gaussian(sigma);
		//scatteringvector = Threevector (scatteringlength_1, scatteringlength_2, scatteringlength_3);
		scatteringvector = Threevector (0., 0., 0.);
		positionnew = positions.back() + scatteringvector + vdrift * scatteringtime * axes.back().normalized();
		axisnew = fTubegeometry->contains(positionnew);
		if (axisnew.compare(control) == true) 
			Nwallcollisions = Nwallcollisions + 1;
		while (axisnew.compare(control) == true){
			
			//scatteringlength_1 = fRandomgenerator->gaussian(sigma); 
			//scatteringlength_2 = fRandomgenerator->gaussian(sigma);
			//scatteringlength_3 = fRandomgenerator->gaussian(sigma);
			//scatteringvector = Threevector (scatteringlength_1, scatteringlength_2, scatteringlength_3);
			scatteringvector = Threevector (0., 0., 0.);
			positionnew = positions.back() + scatteringvector + vdrift * scatteringtime * axes.back().normalized();
			axisnew = fTubegeometry->contains(positionnew);		
		}
		
		times.push_back(tnew);	
		positions.push_back(positionnew);
		axes.push_back(axisnew);
		
		if (fTubegeometry->lastsegmentcontains(positionnew) == true) { 
			wasinlastsegment = true;
			tend = times.back();
		}
	}
	
	if (wasinlastsegment == true && time > tend){
		reachedendoftube = true;
		return positions.back();
	}
	
	else {
		Threevector v;
		double vel;
		int i = Nstart;
		while (times[i] < time)  
			i = i+1;	
		Nstart = i - 1;	
		v = positions[i] + (-1)*positions[i-1];
		vel = v.mag() / (times[i]-times[i-1]);
		pos = positions[i-1] + (time-times[i-1]) * vel * v.normalized();		
		return pos; 
	}
}

void Tubetracking::makeTrack(double t_start, double h){
	if (h > 1e-7) 
		scatteringtime = 1e-7;
	else
		scatteringtime = h;
}

void Tubetracking::reset(){
	Nstart = 0;
}

void Tubetracking::stepDone(double time){
	if (reachedendoftube == true) {
		std::cout << "Position = " << positions.back().toString() << std::endl;
		std::cout << "Axis = " << axes.back().toString() << std::endl;
	}
	
	else {
		int i = Nstart + 1;
		int N = times.size();
		while (times[i] < time)
			i = i + 1;
		if (savetrack == true){
			if(trackparticle.is_open()) {
				for(int j = 0; j < i; j++){
					trackparticle << positions[j][0] << "	" << positions[j][1] << "	" << positions[j][2];
					trackparticle << std::endl;
				}
			}
		}
		for (int j = 0; j < N-i+1; j++){
			times[j] = times[j+i-1];
			axes[j] = axes[j+i-1];
			positions[j] = positions[j+i-1];
		}
		times.resize(N-i+1);
		axes.resize(N-i+1);
		positions.resize(N-i+1);
		Nstart = 0;
	}
}
