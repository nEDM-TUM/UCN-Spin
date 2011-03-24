#include "roots.h"

/**
 * Find intervals in which the function @p f changes sign. This does not find
 * roots without change of sign or two roots that are less than @p dx apart.
 *
 * @param[in]  f     function that is checked for roots
 * @param[in]  start start of the interval in which roots shall be searched
 * @param[in]  stop  end of the interval in which roots shall be searched
 * @param[in]  dx    step size for the finding of roots
 * @param[out] results vector to which the intervals with roots will be added
 */
template <class T>
int Roots::bracketRoots(T &f, double &start, double &stop, double dx, std::vector<std::pair<double,double> > &results) {
	// numer of roots that were found
	int nroots = 0;

	// save last value and position to remember its sign
	double last = 0;
	double last_pos = stop;
	
	// go through the interval
	for (double pos = start; pos <= stop; pos += dx) {
		// get value at current position
		double val = f(pos);

		// check if sign changed
		if (last*val < 0) {
			// add bracket to list of results
			std::pair<double,double> interval(last_pos, pos);
			results.push_back(interval);
			nroots++;
		}

		// save value to compare in next iteration
		last = val;
	}

	return nroots;
}
