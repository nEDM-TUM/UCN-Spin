#include <iostream>
#include "parameters.h"
#include <utility> // make_pair

/**
 * @class Parameters
 * Read and save parameters for the simulation
 *
 * The method #readParameters of this class will read parameters from
 * an input stream of the simple format
 * @verbatim name value @endverbatim
 * and save them for later use.
 *
 * The parameters have to be introduced to the class by calling
 * #expectInt(std::string) or #expectDouble(std::string) for
 * it before calling #readParameters.
 */

Parameters::Parameters()
{

}

Parameters::~Parameters()
{

}

/**
 * Set the value for a parameter.
 * @param name name of the parameter
 * @param value new value of the parameter
 */
void Parameters::add(std::string name, double value)
{
	// don't use insert, it will not replace a default value
	fDoubles[name] = value;
}

/**
 * Set the value for a parameter.
 * @param name name of the parameter
 * @param value new value of the parameter
 */
void Parameters::add(std::string name, int value)
{
	// don't use insert, it will not replace a default value
	fInts[name] = value;
}

/**
 * Tell the parser that we expect a parameter of
 * type int.
 * @param name the name of the expected parameter
 */
void Parameters::expectInt(std::string name) {
	fExpectedInts.insert(name);
}

/**
 * Tell the parser thate we expect a parameter and give
 * default value.
 * @see #expectInt(std::string)
 */
void Parameters::expectInt(std::string name, int def) {
	expectInt(name);
	add(name, def);
}

/**
 * Tell the parser that we expect a parameter of
 * type double.
 * @param name the name of the expected parameter
 */
void Parameters::expectDouble(std::string name) {
	fExpectedDoubles.insert(name);
}

/**
 * Tell the parser thate we expect a parameter and give
 * default value.
 * @see #expectDouble(std::string)
 */
void Parameters::expectDouble(std::string name, double def) {
	expectDouble(name);
	add(name, def);
}

/**
 * Retrieve value of parameter
 */
double Parameters::getDoubleParam(std::string name) const
{
	std::map<std::string,double>::const_iterator it = fDoubles.find(name);
	if(it == fDoubles.end())
	{
		std::cout << "Parameter " << name << " not found, returning 0.0" << std::endl;
		return 0.0;
	}
	return it->second;
}

/**
 * Retrieve value of parameter
 */
int Parameters::getIntParam(std::string name) const
{
	std::map<std::string,int>::const_iterator it = fInts.find(name);
	if(it == fInts.end())
	{
		std::cout << "Parameter " << name << " not found, returning 0" << std::endl;
		return 0;
	}
	return it->second;
}

/**
 * Return total number of parameters in table.
 */
int Parameters::getSize() const
{
	return fDoubles.size() + fInts.size();
}

/**
 * Read parameters from stream.
 * @param in stream from which parameters will be read
 * @see Parameters for description of file format
 */
void Parameters::readParameters(std::istream &in) {
	std::string name;
	double d;
	int i;

	while (!in.eof()) {
		// read name from stream
		in >> name;

		if (name[0] == '#') {
			// ignore comment
		}
		else if (fExpectedInts.count(name) > 0) {
			in >> i; // read integer
			add(name, i); // add to parameters
		}
		else if (fExpectedDoubles.count(name) > 0) {
			in >> d; // read double
			add(name, d); // add to parameters
		}
		else {
			std::cout << "Parameter " << name << " not known!" << std::endl;
		}

		// ignore garbage at end of line
		std::cin.ignore(4192, '\n');
	}
}
