#include "basegeometry.h"
#include "threevector.h"
#include "segment.h"
#include "random.h"
#include "parameters.h"

// Offene Probleme:
// Sicherstellen, dass man beim Auslesen der Datei am Ende der Zeile jeweils angekommen ist?


class TubeGeometry : public Basegeometry {
	public:
// Im Konstruktor wird ein File ausgelesen, das die Informationen über den Schlauch 
// enthält. Der Schlauch soll sich aus Kreisbogenstücken und geraden Stücken
// zusammensetzen. Der Konstruktor erzeugt aus dem File einen Vektor mit Segmenten
// des Schlauchs. Diese Segmente sind entweder Objekte der Klasse Csegment oder Lsegment.
		TubeGeometry(Random*, std:string Tubefile);

// Würfelt einen Startvektor in der Startebene des Schlauchs. Und gibt den ersten
// axialen Vektor zurück. Vorsicht: Diese Funktion verändert die Threevectoren v und x. 
		void initialize(Threevector &v, Threevector &x);
		Threevector boundsCheck(const Threevector &x){
			throw "not implemented";
		};
		void reflect(Threevector &v, const Threevector &x) {
			throw "not implemented";
		};
		double findIntersection(const double t0, const double t1,
				const Polynom &px, const Polynom &py, const Polynom &pz, double eps){
			throw "not implemented";
		};
		bool contains(const Threevector &x, double rootstart) const {
			throw "not implemented";
		};

// 'contains' prüft ob eine Position x im Schlauch enthalten ist oder nicht.
// Wenn ja, gibt sie den axialen Vektor der Projektion der Position auf die Schlauchachse
// zurück. Wenn nicht, dann gibt sie als axialen Vektor den Nullvektor zurück.
// Vorsicht: das zweite Argument 'rootstart' wird während der Funktion verändert.
		Threevector contains(const Threevector &x, double &rootstart);
		bool lastsegmentcontains(const Threevector &x, const double &rootstartsegment);
	
	private:
		Random *rand;
};
