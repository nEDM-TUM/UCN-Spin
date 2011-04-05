all: cylinder roots

cylinder: cylinder.cpp ../objs/cylinder.o
	$(CXX) $(CXXFLAGS) $(LIBS) cylinder.cpp ../objs/cylinder.o -o cylinder

roots: roots.cpp ../roots.h ../objs/polynom.o
	$(CXX) $(CXXFLAGS) $(LIBS) roots.cpp ../objs/polynom.o -o roots

# be careful and let make in .. decide if a prequisite should be rebuilt
../objs/%.o: $(wildcard ../*.cpp) $(wildcard ../*.h)
	make -C .. objs/$*.o

clean:
	rm -f cylinder roots

.PHONY: all clean
