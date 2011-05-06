#include <vector>
#include "random.h"
#include "parameters.h"
#include "basetracking.h"
#include "threevector.h"
#include "tubegeometry.h"

// Offene Probleme:
// Feststellung, ob man am Schlauchende angekommen ist.

// Offene Fragen:
// Wann wird thetubegeometry, Objekt der Klasse Tubegeometry erzeugt?

// Paramter innerhalb dieser Klasse, die noch in die Paramterliste eingefügt werden 
// müssen:
// double mu : in random: exponential(mu) zum Würfeln der Stoßzeit
// double vdrift : Driftgeschwindigkeit der Xe-Atome im Schlauch
// Außerdem muss noch die Berechnung von 'sigma' aus der Streuzeit eingefügt werden.
// Filename, mit den Informationen über den Schlauch.

class Tubetracking : Basetracking {	
	public:
		Tubetracking(Random*);
		
// Ruft Initialize von Tubegeometry auf. Fügt den 1. axialen Vektor zu 'axes' 
// und den Startpunkt zu 'positions' hinzu. Dazu muss 'initialize' einmal zu 
// Beginn aufgerufen werden. 
		virtual void initialize();
		
// Gibt die Position des Xe-Teilchens zu einer bestimmten Zeit 'time' zurück.
		Threevector getPosition(double time);
		virtual void makeTrack(double hmax){
			throw "not implemented";
		};
		void reset();

// Löscht den Inhalt aus den Vektoren 'positions', 'axes', 'times', 'roots' und fügt 
// gleichzeitig die jeweils letzten Werte wieder hinzu.
		void stepDone(double time);
		
		bool reachedendoftube = false;
	
	private:
		Random *rand;
		int Nstart = 0;
		int wasinlastsegment = 0;
// In den Vektoren sollen die bisherigen ausgewürfelten Werte, die zu einer "guten"
// Position geführt haben, gespeichert werden. 
		vector<double> times;
		vector<double> roots;
		vector<Threevector> positions;
		vector<Threevector> axes;
		Tubegeometry *thetubegeomtry;
};
 
