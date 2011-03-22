CXX=g++
CXXFLAGS=-Wall -O3 -march=native -fopenmp
LIBS=-lgsl -lgslcblas -lm
OBJS=main.o bfield.o random.o tracking.o dopr.o derivatives.o parameters.o cylinder.o threevector.o

cylindric : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) -o cylindric

main.o : main.cpp globals.h bfield.h random.h tracking.h dopr.h derivatives.h parameters.h
	$(CXX) $(CXXFLAGS) -c main.cpp

bfield.o : bfield.cpp bfield.h parameters.h
	$(CXX) $(CXXFLAGS) -c bfield.cpp

parameters.o : parameters.cpp parameters.h
	$(CXX) $(CXXFLAGS) -c parameters.cpp

random.o : random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp

tracking.o : tracking.cpp tracking.h random.h bfield.h
	$(CXX) $(CXXFLAGS) -c tracking.cpp

dopr.o : dopr.cpp dopr.h derivatives.h tracking.h
	$(CXX) $(CXXFLAGS) -c dopr.cpp

derivatives.o : derivatives.cpp derivatives.h tracking.h globals.h
	$(CXX) $(CXXFLAGS) -c derivatives.cpp

cylinder.o: cylinder.cpp cylinder.h
	$(CXX) $(CXXFLAGS) -c cylinder.cpp

threevector.o: threevector.cpp threevector.h
	$(CXX) $(CXXFLAGS) -c threevector.cpp

tests/cylinder: tests/cylinder.cpp cylinder.o
	$(CXX) $(CXXFLAGS) $(LIBS) tests/cylinder.cpp cylinder.o -o tests/cylinder

clean : 
	rm -f $(OBJS)
	rm -f cylindric

doc:
	doxygen

.PHONY: clean doc
