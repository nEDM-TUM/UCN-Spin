#include "tubetracking.h"
#include "debug.h"

// Setzt den Zeitstartwert und die Projektion der Startposition auch die Schlauchachse 
// auf 0. 
Tubetracking::Tubetracking(Random *ran, Tubegeometry * geo, Parameters& theParameters)
:Basetracking(ran, geo), rand(ran), fTubegeometry(geo)
{
	times.push_back(0.0);
	//roots.push_back(0.0);
	vdrift = theParameters.getDoubleParam("vdrift");
	mu = theParameters.getDoubleParam("mu");
	sigma = theParameters.getDoubleParam("sigma");
	scatteringtime = 0;
	reachedendoftube = false;
	Nstart = 0;
	wasinlastsegment = 0;
	tend = 0;
	trackparticle.open("track.txt");
}

Tubetracking::~Tubetracking(){
	trackparticle.close();
}

void Tubetracking::initialize(){
        //positions.clear();
        //axes.clear();
        //times.clear();
	Threevector v,x;
	v = Threevector();
	x = Threevector();
	fTubegeometry->initialize(v,x); 
	positions.push_back(x);
	axes.push_back(v);
}
// Prüft zunächst ob die gefragte Zeit schon ausgewürfelt wurde oder noch nicht. 
// Falls nicht, werden so lange "gute" neue Positionen ausgewürfelt, bis die 
// Position zum gefragten Zeitpunkt berechnet werden kann. 
// Zunächst wird eine Streuzeit 'scatteringtime' ausgewürfelt
// 'rootnew' ist die Projektion der neuen Position auf die Schlauchachse. 'rootnew'
// wird von der Funktion 'contains' verändert.
// Zur Ermittlung der neuen Position wird zum gewürfelten Streuvektor 'scatteringvector'
// noch die in der Zeit zurückgelegten Driftlänge in Richtung der vorherigen Achse 
// addiert und mit 'contains' überprüft ob die neue Position innerhalb des Schlauchs 
// liegt oder nicht. Wenn nicht wird für 'axis' der Nullvektor ausgegeben. 
// Falls der neue Ort nicht im Schlauch liegt, wird so lange neu gewürfelt bis er drin ist.
// Zum Schluss werden die neuen Daten in den vier Vektoren gespeichert.
Threevector Tubetracking::getPosition(double time){
	Threevector pos;	
	
	if (time < 1e-50){
		return positions.back();
	}
	
	while (times.back() < time && wasinlastsegment == false) {
		double tnew, scatteringlength_1, scatteringlength_2, scatteringlength_3;
		Threevector positionnew, axisnew, control, scatteringvector;
		//scatteringtime = fRandomgenerator->exponential(mu); 
		tnew = times.back() + scatteringtime;
		//rootnew = roots.back(); 
		control = Threevector (0,0,0);
		
		scatteringlength_1 = scatteringtime * fRandomgenerator->gaussian(sigma); 
		scatteringlength_2 = scatteringtime * fRandomgenerator->gaussian(sigma);
		scatteringlength_3 = scatteringtime * fRandomgenerator->gaussian(sigma);
		scatteringvector = Threevector (scatteringlength_1, scatteringlength_2, scatteringlength_3);
		
		positionnew = positions.back() + scatteringvector + vdrift * scatteringtime * axes.back().normalized();
		//axisnew = fTubegeometry->contains(positionnew, rootnew);
		axisnew = fTubegeometry->contains(positionnew);
		
		while (axisnew.compare(control) == true){
			
			scatteringlength_1 = scatteringtime * fRandomgenerator->gaussian(sigma); 
			scatteringlength_2 = scatteringtime * fRandomgenerator->gaussian(sigma);
			scatteringlength_3 = scatteringtime * fRandomgenerator->gaussian(sigma);
			scatteringvector = Threevector (scatteringlength_1, scatteringlength_2, scatteringlength_3);
			
			positionnew = positions.back() + scatteringvector + vdrift * scatteringtime * axes.back().normalized();
			//axisnew = fTubegeometry->contains(positionnew, rootnew);
			axisnew = fTubegeometry->contains(positionnew);		
		}
		
		times.push_back(tnew);	
		positions.push_back(positionnew);
		axes.push_back(axisnew);
		//roots.push_back(rootnew);
		
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

// Löscht die nicht mehr benötigten Daten aus den vier Vektoren. 
void Tubetracking::stepDone(double time){
	if (reachedendoftube == true) {
		//std::cout << "Position = " << positions.back().toString() << std::endl;
		//std::cout << "Achse = " << axes.back().toString() << std::endl;
		positions.clear();
		axes.clear();
//		roots.clear();
		times.clear();
	}
	else {
		int i = Nstart + 1;
		int N = times.size();
		while (times[i] < time)
			i = i + 1;
		std::cout << "Anzahl der Schritte =" << i << std::endl;
		if (savetrack == true){
			if(trackparticle.is_open()) {
				for (int j = 0; j < i; j++){
					trackparticle << positions[j][0] << "	" << positions[j][1] << "	" << positions[j][2];
					trackparticle << std::endl;
				}
			}
		}
		for (int j = 0; j < N-i+1; j++){
			times[j] = times[j+i-1];
			//roots[j] = roots[j+i-1];
			axes[j] = axes[j+i-1];
			positions[j] = positions[j+i-1];
		}
		times.resize(N-i+1);
		//roots.resize(N-i+1);
		axes.resize(N-i+1);
		positions.resize(N-i+1);
		Nstart = 0;
		std::cout << "Position = " << positions[0].toString() << std::endl;
		//std::cout << "Achse = " << axes[0].toString() << std::endl;
	}
}
