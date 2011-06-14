#ifndef _PARAMETERIZEDFILED_H
#define _PARAMETERIZEDFILED_H

#include "parameters.h"
#include "basetracking.h"
#include "threevector.h"
#include "bfield.h"

class ParameterizedField : public Bfield
{
	public:
		ParameterizedField(const Parameters&, Basetracking* const btr);
		Threevector operator()(const double time) const { return eval(time); };
		Threevector eval(const double time) const;
		
		static const unsigned int ORDERS = 30;

	private:
		void readFieldCoeffs();

		Basetracking* const fTracker;
		double fCoeffs[ORDERS];
};

#endif // _PARAMETERIZEDFILED_H
