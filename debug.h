#ifndef _UCN_DEBUG_H
#define _UCN_DEBUG_H

#include <iostream>

#ifndef NDEBUG
#define debug std::cerr << "(" __FILE__ ":" << __LINE__ << ") : "
#define initialize_debug() std::cerr << std::fixed; std::cerr.precision(15)
#else
#define debug if(false) std::cerr
#define initialize_debug()
#endif

#endif //  _UCN_DEBUG_H
