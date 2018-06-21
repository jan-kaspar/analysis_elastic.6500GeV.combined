HepMC_DIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07
GSL_DIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gsl/1.10
Elegent_DIR = $(shell pwd -P)/../Elegent

INCS = `root-config --cflags` -I$(GSL_DIR)/include -I${Elegent_DIR} -I$(HepMC_DIR)/include
LIBS = `root-config --libs` -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm -L${Elegent_DIR}/lib -lElegent -L$(HepMC_DIR)/lib -lHepMC

all: do_fit printResultTable printResultTableShort

do_fit: do_fit.cc ../HadronicFitModel.h ../classes.h ../command_line_tools.h common.h MethodSimpleFit.h
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -lMinuit -lMinuit2 -Wl,-rpath=${Elegent_DIR}/lib \
		do_fit.cc -o do_fit

printResultTable: printResultTable.cc ../classes.h ../command_line_tools.h
	g++ --std=c++11 `root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` -O3 -g -Wall -Wextra \
		printResultTable.cc -o printResultTable

printResultTableShort: printResultTableShort.cc ../classes.h ../command_line_tools.h
	g++ --std=c++11 `root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` -O3 -g -Wall -Wextra \
		printResultTableShort.cc -o printResultTableShort
