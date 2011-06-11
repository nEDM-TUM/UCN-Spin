#ifndef _EXCEPTIONS_H
#define _EXCEPTIONS_H

class Exception {
	public:
		Exception() : fWhat("unknown exception") {};
		Exception(const char *what) : fWhat(what) {};
		virtual const char *what() const { return fWhat; };
	private:
		const char *fWhat;
};

class EndlessLoop : public Exception {
	public:
		EndlessLoop() : Exception("endless loop") {};
};

#endif // _EXCEPTIONS_H
