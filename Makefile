CXX=g++
CXXFLAGS=-Wall -O3 -march=native -fopenmp
LIBS=-lgsl -lgslcblas -lm
OBJS=main.o bfield.o random.o tracking.o spintracking.o dopr.o derivatives.o parameters.o

cylindric : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) -o cylindric

main.o : main.cpp globals.h bfield.h random.h tracking.h spintracking.h derivatives.h parameters.h
	$(CXX) $(CXXFLAGS) -c main.cpp

bfield.o : bfield.cpp bfield.h parameters.h
	$(CXX) $(CXXFLAGS) -c bfield.cpp

parameters.o : parameters.cpp parameters.h
	$(CXX) $(CXXFLAGS) -c parameters.cpp

random.o : random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp

tracking.o : tracking.cpp tracking.h random.h bfield.h
	$(CXX) $(CXXFLAGS) -c tracking.cpp

spintracking.o : spintracking.cpp spintracking.h dopr.h derivatives.h
	$(CXX) $(CXXFLAGS) -c spintracking.cpp

dopr.o : dopr.cpp dopr.h derivatives.h tracking.h
	$(CXX) $(CXXFLAGS) -c dopr.cpp

derivatives.o : derivatives.cpp derivatives.h tracking.h globals.h
	$(CXX) $(CXXFLAGS) -c derivatives.cpp

clean : 
	rm -f $(OBJS)
	rm -f cylindric

