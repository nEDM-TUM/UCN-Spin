#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>

#include "superpositionfield.h"
#include "basicfields.h"

SuperpositionField::SuperpositionField(Basetracking *t, std::string field_file) : Bfield(t)
{
	readFields(field_file);
}

Threevector SuperpositionField::eval(const double time) const
{
	Threevector ret;

	for (unsigned int i = 0; i < fFields.size(); i++) {
		ret += fFields[i]->eval(time);
	}

	return ret;
}

void SuperpositionField::readFields(std::string fname)
{
	// open file
	std::ifstream f(fname.c_str());

	// check for success
	if (!f.is_open() || !f.good()) {
		std::cerr << "Could not open field file:" << fname << std::endl;
		exit(1);
	}

	// read lines
	while (!f.eof()) {
		std::string sline;
		getline(f, sline); // read one line
		std::stringstream line(sline);

		/// The first token in the line is the type
		std::string type;
		line >> type;

		if (type.size() == 0) {
			/// Empty lines are ignored
		}
		else if (type[0] == '#') {
			/// Comments begin with a # and are only allowed on a
			/// line on their own.
		}
		else if (type == "D") {
			/// D is a magentic dipole. The first three parameters are
			/// the magnetic moment, the next three the position of the
			/// dipole.
			Threevector m;
			Threevector pos;

			line >> m;
			line >> pos;

			fFields.push_back(new DipoleField(fTracker, m, pos));

			std::cout << "Added magnetic dipole with m = " << m << " A m^2" << " and position " << pos << " m" << std::endl;
		}
		else if (type == "E") {
			/// E is a homogenous electric field that causes relativitic effects.
			/// It takes three doubles as arguments
			Threevector E;
			line >> E;

			fFields.push_back(new RelativisticField(fTracker, E));

			std::cout << "Added electric field with E = " << E << " V/m" << std::endl;
		}
		else if (type == "H") {
			/// H is an homogenous magnetic field, just like E. It takes the
			/// vector B as three doubles as argument.
			Threevector B;
			line >> B;

			fFields.push_back(new HomogenousMagneticField(fTracker, B));

			std::cout << "Added homogenous magnetic field B = " << B << " T" << std::endl;
		}
		else if (type == "L") {
			/// L for "laser" is a homogenous magnetic field inside a tube (the laser beam).
			/// It takes the field as three doubles, the height of the center of the tube and
			/// the radius of the tube as parameters
			Threevector B;
			double h, r;
			line >> B;
			line >> h;
			line >> r;

			fFields.push_back(new LaserField(fTracker, B, h, r));

			std::cout << "Added laser field B = " << B << " T" << " inside of cylinder with center at h = " << h
			          << " m and radius of " << r << " m" << std::endl;
		}
		else if (type == "G") {
			/// G is a gradient field and takes only the gradient in z direction
			/// as parameter
			double gradient;
			line >> gradient;

			fFields.push_back(new GradientField(fTracker, gradient));

			std::cout << "Added gradient field B = (r*g/r, r*g/2, g*z) with g = " << gradient << std::endl;
		}
		else {
			/// The program is aborted if an unknown entry is encountered.
			std::cerr << "Can't parse line '" << sline << "' in " << fname << std::endl;
			exit(1);
		}
	}
}

SuperpositionField::~SuperpositionField() {
	for (unsigned int i = 0; i < fFields.size(); i++)
		delete fFields[i];
}
