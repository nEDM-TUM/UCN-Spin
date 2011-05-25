CXX=g++
CXXFLAGS=-Wall -O3 -march=native -fopenmp -ggdb
LIBS=-lgsl -lgslcblas -lm
OBJS=main.o bfield.o random.o dopr.o derivatives.o parameters.o cylinder.o threevector.o basetracking.o equationtracker.o polynom.o lsegment.o csegment.o tubegeometry.o tubetracking.o
TAGFILES=$(shell find . -name "*.cpp" -or -name "*.h")

all: cylindric tags

# This Makefile uses dependency based make.
#
# If you want to compile file.cpp, simply add file.o to OBJS at the top of this
# file.

cylindric : $(OBJS:%.o=objs/%.o)
	$(CXX) $(CXXFLAGS) $(OBJS:%.o=objs/%.o) $(LIBS) -o cylindric

deps/%.d: %.cpp
	$(CXX) $(CXXFLAGS) -MM -MP -MF deps/$*.d -MT objs/$*.o $*.cpp

-include $(OBJS:%.o=deps/%.d)

objs/%.o::
	$(CXX) $(CXXFLAGS) -c $*.cpp -o objs/$*.o

clean : 
	rm -f $(OBJS:%.o=deps/%.d)
	rm -f $(OBJS:%.o=objs/%.o)
	rm -f cylindric
	make -C tests clean

doc:
	doxygen > /dev/null

tags: $(TAGFILES)
	ctags $(TAGFILES)

.PHONY: clean doc all
