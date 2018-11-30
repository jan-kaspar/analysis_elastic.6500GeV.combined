HepMC_DIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07
GSL_DIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gsl/1.10
Elegent_DIR = ../Elegent

INCS = `root-config --cflags` -I$(GSL_DIR)/include -I${Elegent_DIR} -I$(HepMC_DIR)/include
LIBS = `root-config --libs` -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm -L${Elegent_DIR}/lib -lElegent -L$(HepMC_DIR)/lib -lHepMC

all: exploration

exploration: exploration.cc ../HadronicFitModel.h
	g++ --std=c++11 -O3 -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		exploration.cc -o exploration
