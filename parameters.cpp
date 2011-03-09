#include <iostream>
#include "parameters.h"
#include <utility> // make_pair

Parameters::Parameters()
{

}

Parameters::~Parameters()
{

}

void Parameters::add(std::string name, double value)
{
	fDoubles.insert(std::make_pair(name,value));
}

void Parameters::add(std::string name, int value)
{
	fInts.insert(std::make_pair(name,value));
}

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

int Parameters::getSize() const
{
	return fDoubles.size() + fInts.size();
}
