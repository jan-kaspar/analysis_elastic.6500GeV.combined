HepMC_DIR = /afs/cern.ch/cms/slc6_amd64_gcc481/external/hepmc/2.06.07
GSL_DIR = /afs/cern.ch/cms/slc6_amd64_gcc481/external/gsl/1.10/
Elegent_DIR = ../Elegent

INCS = `root-config --cflags` -I$(GSL_DIR)/include -I${Elegent_DIR} -I$(HepMC_DIR)/include
LIBS = `root-config --libs` -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm -L${Elegent_DIR}/lib -lElegent -L$(HepMC_DIR)/lib -lHepMC

all: simulation normalize_and_fit process_fits

simulation: simulation.cc ../HadronicFitModel.h
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		simulation.cc -o simulation

normalize_and_fit: normalize_and_fit.cc ../HadronicFitModel.h fit.h
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -lMinuit -lMinuit2 -Wl,-rpath=${Elegent_DIR}/lib \
		normalize_and_fit.cc -o normalize_and_fit

process_fits: process_fits.cc
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		process_fits.cc -o process_fits
