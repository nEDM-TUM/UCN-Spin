all: cylinder roots

cylinder: cylinder.cpp
	make -C .. cylinder.o
	$(CXX) $(CXXFLAGS) $(LIBS) cylinder.cpp ../cylinder.o -o cylinder

roots: roots.cpp ../roots.h
	make -C .. polynom.o
	$(CXX) $(CXXFLAGS) $(LIBS) roots.cpp ../polynom.o -o roots

clean:
	rm -f cylinder

.PHONY: all clean
