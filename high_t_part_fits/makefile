INCS = `root-config --cflags`
LIBS = `root-config --libs`

all: do_fit

do_fit: do_fit.cc
	g++ --std=c++11 -O3 -g -Wall -Wextra $(INCS) $(LIBS) -lMinuit -lMinuit2 \
		do_fit.cc -o do_fit
