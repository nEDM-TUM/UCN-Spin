LIBS=-lgsl -lgslcblas -lm
CXXFLAGS=-ggdb

all: cylinder roots polynom tube

cylinder: cylinder.cpp ../objs/cylinder.o ../objs/threevector.o ../objs/random.o ../objs/polynom.o ../objs/basetracking.o ../objs/equationtracker.o ../interpolationpolynomial.h
	$(CXX) $(CXXFLAGS) $(LIBS) cylinder.cpp ../objs/cylinder.o ../objs/random.o  ../objs/polynom.o ../objs/basetracking.o ../objs/equationtracker.o ../objs/threevector.o -o cylinder

roots: roots.cpp ../roots.h ../objs/polynom.o
	$(CXX) $(CXXFLAGS) $(LIBS) roots.cpp ../objs/polynom.o -o roots

polynom: polynom.cpp ../objs/polynom.o
	$(CXX) $(CXXFLAGS) $(LIBS) polynom.cpp ../objs/polynom.o -o polynom

tube: tube.cpp ../objs/tubegeometry.o ../objs/threevector.o 
	$(CXX) $(CXXFLAGS) $(LIBS) tube.cpp ../objs/tubegeometry.o ../objs/threevector.o -o tube

# be careful and let make in .. decide if a prequisite should be rebuilt
../objs/%.o: $(wildcard ../*.cpp) $(wildcard ../*.h)
	make -C .. objs/$*.o

clean:
	rm -f cylinder roots

.PHONY: all clean
