HepMC_DIR = /afs/cern.ch/cms/slc6_amd64_gcc481/external/hepmc/2.06.07
GSL_DIR = /afs/cern.ch/cms/slc6_amd64_gcc481/external/gsl/1.10/
Elegent_DIR = $(shell pwd -P)/../Elegent

INCS = `root-config --cflags` -I$(GSL_DIR)/include -I${Elegent_DIR} -I$(HepMC_DIR)/include
LIBS = `root-config --libs` -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm -L${Elegent_DIR}/lib -lElegent -L$(HepMC_DIR)/lib -lHepMC

all: impactParameterDistributions

impactParameterDistributions: impactParameterDistributions.cc NumericalIntegration.h
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -lMinuit -lMinuit2 -Wl,-rpath=${Elegent_DIR}/lib \
		impactParameterDistributions.cc -o impactParameterDistributions
