/**
 * @file
 * This file provides debugging facilities for the code. You can use debug as output
 * stream like cout and cerr, but it will be deactivated, if NDEBUG or NO_DEBUG_PRINTS
 * is set. It will also add the filename and linenumber of the code that generated the
 * debug print to the output. Don't forget to add std::endl at the end of any debug
 * statement.
 */

#ifndef _UCN_DEBUG_H
#define _UCN_DEBUG_H

#include <iostream>

/// Debug prints are automatically deactivated if NDEBUG (which deactivates
/// assert) is set.
#ifdef NDEBUG
#define NO_DEBUG_PRINTS
#endif

/// By setting NO_DEBUG_PRINTS, you can deactivate debug prints without
/// deactivating asserts.
#ifndef NO_DEBUG_PRINTS
#define debug std::cerr << "(" __FILE__ ":" << __LINE__ << ") : "
#define initialize_debug() std::cerr << std::fixed; std::cerr.precision(15)
#else
#define debug if(false) std::cerr
#define initialize_debug()
#endif

#endif //  _UCN_DEBUG_H
