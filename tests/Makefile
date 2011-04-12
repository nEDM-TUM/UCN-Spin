LIBS=-lgsl -lgslcblas -lm

all: cylinder roots polynom

cylinder: cylinder.cpp ../objs/cylinder.o ../objs/threevector.o ../objs/random.o ../objs/polynom.o
	$(CXX) $(CXXFLAGS) $(LIBS) cylinder.cpp ../objs/cylinder.o ../objs/random.o  ../objs/polynom.o ../objs/threevector.o -o cylinder

roots: roots.cpp ../roots.h ../objs/polynom.o
	$(CXX) $(CXXFLAGS) $(LIBS) roots.cpp ../objs/polynom.o -o roots

polynom: polynom.cpp ../objs/polynom.o
	$(CXX) $(CXXFLAGS) $(LIBS) polynom.cpp ../objs/polynom.o -o polynom

# be careful and let make in .. decide if a prequisite should be rebuilt
../objs/%.o: $(wildcard ../*.cpp) $(wildcard ../*.h)
	make -C .. objs/$*.o

clean:
	rm -f cylinder roots

.PHONY: all clean
