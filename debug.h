#ifndef _UCN_DEBUG_H
#define _UCN_DEBUG_H

#include <iostream>

#ifndef NDEBUG
#define debug std::cerr << "(" __FILE__ ":" << __LINE__ << ") : "
#else
#define debug if(false) std::cerr
#endif

#endif //  _UCN_DEBUG_H
