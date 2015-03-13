# LBImproved C++ Library

This library comes in the form of one short C++ header file. The documentation
is in the C++ comments and in this file.


# Key features

1) Fast Dynamic Time Warping nearest neighbor retrieval.

2) Persistence

3) External-memory: you need only a constant amount of RAM


# PREREQUISITES 

1) You must first build and install the spatial index library 
(http://research.att.com/~marioh/spatialindex/index.html)
I built this software with release 1.3.2 - May 23rd, 2008. 

2) While not strictly necessary, SWIG (http://www.swig.org) is strongly recommended.
I interact with the library using swig and python.

3) If you are using SWIG, Python is recommended. I have used Python 2.5.

# OPERATING SYSTEM 

I built and ran this software with Mac OS 10.4. It also builds under Linux
if you have Python 2.5 and swig installed. It should be possible to
use  any other Unix-like operating system, or even Windows. 

# BUILD 

type "make"


# TESTING 

type "python unitesting.py"


# USAGE 


    import dtw
    constraint = 0.1
    n = 128
    c =int(constraint*n)
    rtree=dtw.TimeSeriesTree("mytmpfile.bin",c,reducdim)
    # randomwalk(n) return a size n array
    for i in xrange(1000):
       rtree.add(randomwalk(n))
    x = randomwalk(n)
    for mode in [rtree.LINEAR, rtree.TREE]:
     for algo in [ rtree.NAIVE,rtree.LB_KEOGH, rtree.LB_IMPROVED]:
       rtree.getNearestNeighborCost(x,algo,mode) 
    rtree.close()
    # to reopen the tree, just do this:
    rtree=dtw.TimeSeriesTree("mytmpfile.bin")
