all: cylinder

cylinder: cylinder.cpp
	make -C .. cylinder.o
	$(CXX) $(CXXFLAGS) $(LIBS) cylinder.cpp ../cylinder.o -o cylinder

clean:
	rm -f cylinder

.PHONY: all clean
