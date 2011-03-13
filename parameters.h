#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <map>
#include <set>
#include <string>
#include <iostream>

class Parameters
{
	public:
	
	Parameters();
	~Parameters();
	void add(std::string name, double value);
	void add(std::string name, int value);
	void expectInt(std::string name);
	void expectInt(std::string name, int def);
	void expectDouble(std::string name);
	void expectDouble(std::string name, double def);
	double getDoubleParam(std::string name) const;
	int getIntParam(std::string name) const;
	int getSize() const;
	void readParameters(std::istream &in);

	private:
	
	std::map<std::string,double> fDoubles;
	std::map<std::string,int> fInts;
	std::set<std::string> fExpectedDoubles;
	std::set<std::string> fExpectedInts;
};
#endif
