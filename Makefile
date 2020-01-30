all: unittesting benchmark example
OSNAME = $(shell uname | tr "[:upper:]" "[:lower:]")
SHAREDNAME=$(shell if [  $(OSNAME) = "darwin" ]; then echo -n "   -bundle -flat_namespace -undefined suppress"; else echo -n "-shared";fi )


benchmark: benchmarks/benchmark.cpp include/dtw.h 
	$(CXX) -O2 -Wall -Wold-style-cast  -Woverloaded-virtual -o benchmark benchmarks/benchmark.cpp  -Iinclude

example: examples/example.cpp include/dtw.h 
	$(CXX) -O2 -Wall -Wold-style-cast  -Woverloaded-virtual -o example examples/example.cpp  -Iinclude


unittesting: tests/unittesting.cpp include/dtw.h 
	$(CXX) -g3 -Wall -Wold-style-cast  -Woverloaded-virtual -o unittesting tests/unittesting.cpp -Iinclude


clean :
	rm -f benchmark unittesting example
