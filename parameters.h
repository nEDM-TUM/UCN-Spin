#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <map>
#include <string>

class Parameters
{
	public:
	
	Parameters();
	~Parameters();
	void add(std::string name, double value);
	void add(std::string name, int value);
	double getDoubleParam(std::string name) const;
	int getIntParam(std::string name) const;
	int getSize() const;

	private:
	
	std::map<std::string,double> fDoubles;
	std::map<std::string,int> fInts;
};
#endif
