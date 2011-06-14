#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "parameterizedfield.h"

ParameterizedField::ParameterizedField(const Parameters&, Basetracking* const btr) :
	fTracker(btr)
{
	readFieldCoeffs();
}

Threevector ParameterizedField::eval(double time) const
{
	const Threevector position = fTracker->getPosition(time);
	const double x = position[0];
	const double y = position[1];
	const double z = position[2];
	const double s = x*x + y*y + z*z;

	// basis function values
	double bfunc[ORDERS][3];

	// calclulate result
	bfunc[0][0]=0.;
	bfunc[0][1]=1.;
	bfunc[0][2]=0.;

	bfunc[1][0]=0.;
	bfunc[1][1]=0.;
	bfunc[1][2]=1.;

	bfunc[2][0]=1.;
	bfunc[2][1]=0.;
	bfunc[2][2]=0.;

	bfunc[3][0]=-2.*x;
	bfunc[3][1]=4.*y;
	bfunc[3][2]=-2.*z;

	bfunc[4][0]=0.;
	bfunc[4][1]=z;
	bfunc[4][2]=y;

	bfunc[5][0]=-2.*x;
	bfunc[5][1]=-2.*y;
	bfunc[5][2]=4.*z;

	bfunc[6][0]=y;
	bfunc[6][1]=x;
	bfunc[6][2]=0.;

	bfunc[7][0]=z;
	bfunc[7][1]=0.;
	bfunc[7][2]=x;

	bfunc[8][0]=-6.*x*y;
	bfunc[8][1]=-3.*s+9.*pow(y,2);
	bfunc[8][2]=-6.*y*z;

	bfunc[9][0]=-2.*x*z;
	bfunc[9][1]=8.*y*z;
	bfunc[9][2]=-pow(x,2)+4.*pow(y,2)-3.*pow(z,2);

	bfunc[10][0]=-2.*x*y;
	bfunc[10][1]=-pow(x,2)-3.*pow(y,2)+4.*pow(z,2);
	bfunc[10][2]=8.*y*z;

	bfunc[11][0]=-6.*x*z;
	bfunc[11][1]=-6.*y*z;
	bfunc[11][2]=-3.*s+9.*pow(z,2);

	bfunc[12][0]=-3.*pow(x,2)+4.*pow(y,2)-pow(z,2);
	bfunc[12][1]=8.*x*y;
	bfunc[12][2]=-2.*x*z;

	bfunc[13][0]=y*z;
	bfunc[13][1]=x*z;
	bfunc[13][2]=x*y;

	bfunc[14][0]=-3.*pow(x,2)-pow(y,2)+4.*pow(z,2);
	bfunc[14][1]=-2.*x*y;
	bfunc[14][2]=8.*x*z;

	bfunc[15][0]=12.*s*x-60.*x*pow(y,2);
	bfunc[15][1]=-48.*s*y+80.*pow(y,3);
	bfunc[15][2]=12.*s*z-60.*pow(y,2)*z;

	bfunc[16][0]=-6.*x*y*z;
	bfunc[16][1]=-3.*s*z+15.*pow(y,2)*z;
	bfunc[16][2]=-3.*s*y+7.*pow(y,3)-6.*y*pow(z,2);

	bfunc[17][0]=4.*s*x-10.*x*pow(y,2)-10.*x*pow(z,2);
	bfunc[17][1]=-6.*s*y-10.*pow(y,3)+60.*y*pow(z,2);
	bfunc[17][2]=-6.*s*z+60.*pow(y,2)*z-10.*pow(z,3);

	bfunc[18][0]=-6.*x*y*z;
	bfunc[18][1]=-3.*s*z-6.*pow(y,2)*z+7.*pow(z,3);
	bfunc[18][2]=-3.*s*y+15.*y*pow(z,2);

	bfunc[19][0]=12.*s*x-60.*x*pow(z,2);
	bfunc[19][1]=12.*s*y-60.*y*pow(z,2);
	bfunc[19][2]=-48.*s*z+80.*pow(z,3);

	bfunc[20][0]=-3.*s*y-6.*pow(x,2)*y+7.*pow(y,3);
	bfunc[20][1]=-3.*s*x+15.*x*pow(y,2);
	bfunc[20][2]=-6.*x*y*z;

	bfunc[21][0]=-(s*z)-2.*pow(x,2)*z+7.*pow(y,2)*z;
	bfunc[21][1]=12.*x*y*z;
	bfunc[21][2]=-(s*x)+7.*x*pow(y,2)-2.*x*pow(z,2);

	bfunc[22][0]=-(s*y)-2.*pow(x,2)*y+7.*y*pow(z,2);
	bfunc[22][1]=-(s*x)-2.*x*pow(y,2)+7.*x*pow(z,2);
	bfunc[22][2]=12.*x*y*z;

	bfunc[23][0]=-3.*s*z-6.*pow(x,2)*z+7.*pow(z,3);
	bfunc[23][1]=-6.*x*y*z;
	bfunc[23][2]=-3.*s*x+15.*x*pow(z,2);

	bfunc[24][0]=60.*s*x*y-140.*x*pow(y,3);
	bfunc[24][1]=15.*pow(s,2)-150.*s*pow(y,2)+175.*pow(y,4);
	bfunc[24][2]=60.*s*y*z-140.*pow(y,3)*z;

	bfunc[25][0]=4.*s*x*z-28.*x*pow(y,2)*z;
	bfunc[25][1]=-24.*s*y*z+56.*pow(y,3)*z;
	bfunc[25][2]=pow(s,2)-14.*s*pow(y,2)+21.*pow(y,4)+4.*s*pow(z,2)-28.*pow(y,2)*pow(z,2);

	bfunc[26][0]=12.*s*x*y-14.*x*pow(y,3)-42.*x*y*pow(z,2);
	bfunc[26][1]=3.*pow(s,2)-9.*s*pow(y,2)-14.*pow(y,4)-21.*s*pow(z,2)+147.*pow(y,2)*pow(z,2);
	bfunc[26][2]=-30.*s*y*z+112.*pow(y,3)*z-42.*y*pow(z,3);

	bfunc[27][0]=12.*s*x*z-42.*x*pow(y,2)*z-14.*x*pow(z,3);
	bfunc[27][1]=-30.*s*y*z-42.*pow(y,3)*z+112.*y*pow(z,3);
	bfunc[27][2]=3.*pow(s,2)-21.*s*pow(y,2)-9.*s*pow(z,2)+147.*pow(y,2)*pow(z,2)-14.*pow(z,4);

	bfunc[28][0]=4.*s*x*y-28.*x*y*pow(z,2);
	bfunc[28][1]=pow(s,2)+4.*s*pow(y,2)-14.*s*pow(z,2)-28.*pow(y,2)*pow(z,2)+21.*pow(z,4);
	bfunc[28][2]=-24.*s*y*z+56.*y*pow(z,3);

	bfunc[29][0]=60.*s*x*z-140.*x*pow(z,3);
	bfunc[29][1]=60.*s*y*z-140.*y*pow(z,3);
	bfunc[29][2]=15.*pow(s,2)-150.*s*pow(z,2)+175.*pow(z,4);

	Threevector ret(0, 0, 0);

	// add orders and output as threevector
	for (unsigned int o = 0; o < ORDERS; o++) {
		for (unsigned int i = 0; i < 3; i++)
			ret[i] += fCoeffs[o]*bfunc[o][i];
	}

	return ret;
}

void ParameterizedField::readFieldCoeffs() {
	std::ifstream coeff_file("field_coeffs.dat");

	unsigned int i = 0; // number of coeff that is read

	if (!coeff_file.is_open() || !coeff_file.good()) {
		std::cerr << "Could not open coefficient file!" << std::endl;
		exit(1);
	}

	// read coeffs from file
	while (!coeff_file.eof() && i < ORDERS) {
		coeff_file >> fCoeffs[i];
		i++;
	}

	// fill rest with zeros
	while (i < ORDERS) {
		fCoeffs[i] = 0;
		i++;
	}

	// close input file
	coeff_file.close();
}
