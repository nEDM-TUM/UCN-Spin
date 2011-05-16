#include "tubetracking.h"

// Setzt den Zeitstartwert und die Projektion der Startposition auch die Schlauchachse 
// auf 0. 
Tubetracking::Tubetracking(Random *ran)
:rand(ran){
	std::string geometryfilename = "thetube.txt";
	thetubegeometry = new Tubegeometry(ran, geometryfilename);
	times.push_back(double 0);
	roots.push_back(double 0);
	vdrift = theParameters.getDoubleParam("vdrift");
	mu = theParameters.getDoubleParam("mu");
	sigma = theParameters.getDoubleParam("sigma");
}

void Tubetracking::initialize(){
	Threevector v,x;
	v = Threevector();
	x = Threevector();
	thetubegeometry->initialize(v,x); 
	positions.push_back(x);
	axis.push_back(v);
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
	while (times.back() < time && t_end != 0) {
		double tnew, rootnew, simga, scatteringtime, scatteringlength_1, scatteringlength_2, scatteringlength_3;
		Threevector positionnew, axisnew, control, scatteringvector;
		scatteringtime = random->exponential(mu); 
		tnew = times.back() + scatteringtime;
		rootnew = 0; 
		control = Threevector ();
		scatteringlength_1 = scatteringtime * random->gaussian(sigma); 
		scatteringlength_2 = scatteringtime * random->gaussian(sigma);
		scatteringlength_3 = scatteringtime * random->gaussian(sigma);
		scatteringvector = Threevector (scatteringlength_1, scatteringlength_2, scatteringlength_3);
		positionnew = positions.back() + scatteringvector + vdrift * scatteringtime * axes.back().normalize();
		axisnew = thetubegeometry->contains(positionnew, roots.back());
		while (axisnew == control){
			if (thetubegeometry->lastsegmentcontains(positionnew) == true){
				reachedendoftube1 = true;
				axisnew = 
			}
			else {
				scatteringlength_1 = scatteringtime * random->gaussian(sigma); 
				scatteringlength_2 = scatteringtime * random->gaussian(sigma);
				scatteringlength_3 = scatteringtime * random->gaussian(sigma);
				scatteringvector = Threevector (scatteringlength_1, scatteringlength_2, scatteringlength_3);
				positionnew = positions.back() + scatteringvector + vdrift * scatteringtime * axes.back().normalize();
				axisnew = thetubegeomtry->contains(positionnew, roots.back());
			} 
		times.push_back(tnew);	
		positions.push_back(positionnew);
		axes.push_back(axisnew);
		roots.push_back(rootnew);
		if (thetubegeomtry->lastsegmentcontains(positionnew) == true) {
			t_end = times.back();
		}
	}
	
// Ist die gefragte Zeit dann ausgewürfelt, dann wird zuerst der geeignete Schritt
// ermittelt und über lineare Interpolation zwischen den zwei umliegenden Punkten der
// gefragte Ort ermittelt.
	if (t_end != 0 && time > t_end){
		wasinlastsegment = true;
		return positions.back()
	}
	else {
		double vel;
		int i = Nstart;
		while (times[i] < time)  
			i = i+1;	
		Nstart = i - 1;	
		v = positions[i]-positions[i-1];
		vel = v.mag() / (times[i]-times[i-1]);
		pos = positions[i-1] + (time-times[i-1]) * vel * v.normalize();		
		return pos; 
}

void Tubetracking::reset(){
	Nstart = 0;
}

// Löscht die nicht mehr benötigten Daten aus den vier Vektoren. 
void Tubetracking::stepDone(double time){
	int i = Nstart
	int N = times.size()
	while (times[i] < time)  
		i = i + 1;
	for (int j = 0, j < N-i+1, j++){
		times[j] = times[j+i-1];
		roots[j] = roots[j+i-1];
		axes[j] = axes[j+i-1];
		positions[j] = positions[j+i-1];
	}
	times.resize(N-i+1);
	roots.resize(N-i+1);
	axes.resize(N-i+1);
	positions.resize(N-i+1);
	Nstart = 0;
}
