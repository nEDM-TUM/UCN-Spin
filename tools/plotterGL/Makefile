LDFLAGS=-lglut `pkg-config --libs glu glew`
CXXFLAGS=-ggdb -Wall

plotter: plotter.o

plotter.o: plotter.cpp plotter.h

clean:
	rm plotter.o plotter

.PHONY: clean
