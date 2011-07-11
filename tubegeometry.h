#ifndef TUBEGEOMETRY_H
#define TUBEGEOMETRY_H
#include "basegeometry.h"
#include "threevector.h"
#include "lsegment.h"
#include "csegment.h"
#include "random.h"
#include "parameters.h"
#include <string>
#include <vector>


class Tubegeometry : public Basegeometry {
	public:
		Tubegeometry(Random*, std::string Tubefile, std::string Tubefilemathematica, Parameters& theParameters);
// Vorsicht: Diese Funktion verändert die Threevectoren v und x. 
		void initialize(Threevector &v, Threevector &x);
		bool boundsCheck(const Threevector &x){
			throw "not implemented";
		};
		void reflect(Threevector &v, const Threevector &x) {
			throw "not implemented";
		};
		double findIntersection(const double t0, const double t1,
				const Polynom &px, const Polynom &py, const Polynom &pz, double eps){
			throw "not implemented";
		};
		bool contains(const Threevector &x) const {
			throw "not implemented";
		};
		virtual Threevector contains(const Threevector &x);
		virtual bool lastsegmentcontains(const Threevector &x);
	
	private:
		Random *random;
		double radiustube;
		std::vector<Segment*> Segments;
		std::vector<Lsegment> Lsegments;
		std::vector<Csegment> Csegments;
};
#endif
