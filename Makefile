CC=g++
CFLAGS=-c -Wall -std=c++11 -O3 -g
INCLUDES=-I. -I./analysis -I$(shell root-config --incdir)
#(root-config --incdir)
LIBDIRS=-L/usr/lib
LIBS=
LDFLAGS=$(shell root-config --glibs)

SRC=TrackModelPhysical.cxx FitLayer.cxx TrackletConstructor.cxx Tracker.cxx reconstruction.cxx Geo.cxx Cuts.cxx analysis/PolynomFit.cxx DataStructures.cxx

SOURCES=$(SRC) analysis/TrackletFitTest.cxx TrackletProlongation.cxx analysis/TrackletConstructorMC.cxx TrackletCleaner.cxx \
TrackletEfficiency.cxx ReadEvent.cxx TrackletConstructorPrim.cxx analysis/TrackletConstructorPrimMC.cxx \
analysis/AnalyseGeometry.cxx

HEADERS=$(SRC:.cxx=.h) util.h
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=reco

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

%.o: %.cxx $(HEADERS) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ 


clean:
	rm reco $(OBJECTS)

cleanall:
	rm reco $(OBJECTS) *.~ *.root
