all: printResultTable

printResultTable: printResultTable.cc ../classes.h ../../command_line_tools.h
	g++ --std=c++11 `root-config --libs` -lMinuit -lMinuit2 `root-config --cflags` -O3 -g -Wall -Wextra \
		printResultTable.cc -o printResultTable
