#ifndef _SUPERPOSITIONFIELD_H
#define _SUPERPOSITIONFIELD_H

#include <string>
#include <vector>

#include "bfield.h"
#include "basetracking.h"

class SuperpositionField : public Bfield
{
	public:
		SuperpositionField(Basetracking *t, std::string field_file);
		virtual ~SuperpositionField();
		Threevector eval(const double time) const;
	private:
		void readFields(std::string filename);

		std::vector<Bfield*> fFields;
};

#endif // _SUPERPOSITIONFIELD_H
