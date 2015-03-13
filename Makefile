all: swig
OSNAME = $(shell uname | tr "[:upper:]" "[:lower:]")
SHAREDNAME=$(shell if [  $(OSNAME) = "darwin" ]; then echo -n "   -bundle -flat_namespace -undefined suppress"; else echo -n "-shared";fi )
package:
	zip -9 lbimproved_`date +%Y-%m-%d`.zip dtw.h Makefile rtreebased.h querystrategy.h unitesting.py timeseries.i unittesting.cpp README data/*py data/READ*

swig:
	# this assumes you have swig on your system
	swig -Wall -c++ -python timeseries.i
	# your python include directory may vary
	g++ -c -fPIC -O2 timeseries_wrap.cxx   -I/sw/include -I/sw/include/python2.5 -I/usr/include/python2.5/ -I/usr/local/include/spatialindex
	g++  $(SHAREDNAME)  timeseries_wrap.o /usr/local/lib/libspatialindex.a -o _dtw.so 

unittesting: unittesting.cpp dtw.h rtreebased.h
	g++ -g3 -Wall -Wold-style-cast  -Woverloaded-virtual -o unittesting unittesting.cpp
