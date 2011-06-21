#ifndef _TIMEOUT_H
#define _TIMEOUT_H

#include <ctime>
#include "exceptions.h"

class TimeoutException : public Exception {
	public:
		TimeoutException() : Exception("time out") {};
};

class Timeout {
	public:
		explicit Timeout(time_t seconds) : fTimeout(seconds) { reset(); };
		void reset() { fStartTime = time(NULL); };
		void check()
		{
			if (time(NULL) - fStartTime > fTimeout)
				throw TimeoutException();
		};
	private:
		time_t fTimeout;
		time_t fStartTime;
};

#endif // _TIMEOUT_H
