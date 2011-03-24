all: cylinder roots

cylinder: cylinder.cpp
	make -C .. cylinder.o
	$(CXX) $(CXXFLAGS) $(LIBS) cylinder.cpp ../cylinder.o -o cylinder

roots: roots.cpp
	make -C .. roots.o
	$(CXX) $(CXXFLAGS) $(LIBS) roots.cpp ../roots.o -o roots

clean:
	rm -f cylinder

.PHONY: all clean
