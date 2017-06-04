# LBImproved C++ Library
[![Build Status](https://travis-ci.org/lemire/lbimproved.png)](https://travis-ci.org/lemire/lbimproved)

This library comes in the form of one short C++ header file. The documentation
is in the C++ comments and in this file.


# Key feature

1) Fast Dynamic Time Warping nearest neighbor retrieval.

2) Implementations of LB Koegh and LB Improved

3) Companion to the following paper :

Daniel Lemire, Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound, Pattern Recognition 42 (9), pages 2169-2180, 2009. 
http://arxiv.org/abs/0811.3301

Comments about this paper by Keogh's team: 

     To our knowledge, there is only one paper that
     offers a plausible speedup based on a tighter 
     lower bound—Lemire (2009) suggests a mean speedup 
     of about 1.4 based on a tighter bound. 
     These results are reproducible, and testing on 
     more general data sets we obtained similar 
     results (...) (Wang et al. 2013)


# BUILD 

type "make"

    make
    ./unittesting
    ./benchmark
    ./example

# Simple code example

See ``example.cpp``.

# Other libraries
 
*  [dtwclust](https://github.com/asardaes/dtwclust) is an  R Package for Time Series Clustering Along with Optimizations for DTW
