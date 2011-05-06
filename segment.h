# include "threevector.h"
#include "parameters.h"

class Segment {
	public:
		virtual Threevector getposition (double tau) = 0;
		virtual Threevector segmentcontains(const Threevector &x, double &rootstartsegment);
		virtual Threevector axis(double tau) = 0;
		virtual Threevector start() = 0;
};

// Kreissegmente:
// stattet die Objekte mit den erforderlichen Größen aus, um das Segment paramtetrisieren
// zu können, und um festzustellen ob das Segment einen Vektor x enthält oder nicht.
class Csegment : public Segments {
	Threevector a;
	Threevector b;
	Threevector centre;
	Threevector start;
	double radius;
	double t_max;
	Csegment (Threevector s, Threevector v, Threevector n, double r, double t);
	Threevector getposition (double tau);
	Threevector axis(double tau);
	Treevector startpoint();
	double derivdist(const Threevector &x, double tau);
	double secderivdist(const Threevector &x, double tau);
	double rootNewton(const Threevector &x, double rootstartsegment);
// Achtung 'rootstartsegment' wird hier verändert!!
	Threevector segmentcontains(const Threevector &x, double &rootstartsegment); 
};

// Das Gleiche mit geraden Segmenten.
class Lsegment : public Segments {
	Threevector	start;
	Threevector direction;
	double t_max;
	Lsegment (Threevector s, Threevector v, double t);
	Threevector getposition (double tau);
	Threevector axis();
	Threevector axis(double tau);
	Threevector startpoint();
// Achtung 'rootstartsegmet' wird her verändert!!!
	Threevector segmentcontains(const Threevector &x, double &rootstartsegment);
};
