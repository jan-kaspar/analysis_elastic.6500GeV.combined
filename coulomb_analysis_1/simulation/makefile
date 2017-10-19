HepMC_DIR = /afs/cern.ch/cms/slc6_amd64_gcc481/external/hepmc/2.06.07
GSL_DIR = /afs/cern.ch/cms/slc6_amd64_gcc481/external/gsl/1.10/
Elegent_DIR = ../../Elegent

INCS = `root-config --cflags` -I$(GSL_DIR)/include -I${Elegent_DIR} -I$(HepMC_DIR)/include
LIBS = `root-config --libs` -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm -L${Elegent_DIR}/lib -lElegent -L$(HepMC_DIR)/lib -lHepMC

all: build_models build_uncertainties simulation

build_models: build_models.cc ../../HadronicFitModel.h
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		build_models.cc -o build_models

build_uncertainties: build_uncertainties.cc
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		build_uncertainties.cc -o build_uncertainties

simulation: simulation.cc ../command_line_tools.h
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -Wl,-rpath=${Elegent_DIR}/lib \
		simulation.cc -o simulation
