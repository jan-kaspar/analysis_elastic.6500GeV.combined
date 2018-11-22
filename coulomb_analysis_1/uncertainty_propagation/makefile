HepMC_DIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07
GSL_DIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gsl/1.10
Elegent_DIR = $(shell pwd -P)/../../Elegent

INCS = `root-config --cflags` -I$(GSL_DIR)/include -I${Elegent_DIR} -I$(HepMC_DIR)/include
LIBS = `root-config --libs` -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm -L${Elegent_DIR}/lib -lElegent -L$(HepMC_DIR)/lib -lHepMC

all: build_uncertainties modification process_fits

build_uncertainties: build_uncertainties.cc
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		build_uncertainties.cc -o build_uncertainties

modification: modification.cc
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		modification.cc -o modification

process_fits: process_fits.cc
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		process_fits.cc -o process_fits
