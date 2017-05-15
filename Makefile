all: unittesting benchmark example
OSNAME = $(shell uname | tr "[:upper:]" "[:lower:]")
SHAREDNAME=$(shell if [  $(OSNAME) = "darwin" ]; then echo -n "   -bundle -flat_namespace -undefined suppress"; else echo -n "-shared";fi )
package:
	zip -9 lbimproved_`date +%Y-%m-%d`.zip dtw.h Makefile rtreebased.h querystrategy.h unitesting.py timeseries.i unittesting.cpp README data/*py data/READ*

benchmark: benchmark.cpp dtw.h 
	$(CXX) -O2 -Wall -Wold-style-cast  -Woverloaded-virtual -o benchmark benchmark.cpp

example: example.cpp dtw.h 
	$(CXX) -O2 -Wall -Wold-style-cast  -Woverloaded-virtual -o example example.cpp


unittesting: unittesting.cpp dtw.h 
	$(CXX) -g3 -Wall -Wold-style-cast  -Woverloaded-virtual -o unittesting unittesting.cpp


clean :
	rm -f benchmark unittesting example
