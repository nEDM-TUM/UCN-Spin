#include <vector>
#include "random.h"
#include "parameters.h"
#include "basetracking.h"
#include "threevector.h"
#include "tubegeometry.h"

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
		double vdrift, mu, sigma;
		Random *rand;
		int Nstart = 0;
		int wasinlastsegment = 0;
		int t_end = 0;
// In den Vektoren sollen die bisherigen ausgewürfelten Werte, die zu einer "guten"
// Position geführt haben, gespeichert werden. 
		vector<double> times;
		vector<double> roots;
		vector<Threevector> positions;
		vector<Threevector> axes;
		Tubegeometry *thetubegeomtry;
};
 
