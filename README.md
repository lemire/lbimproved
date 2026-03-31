# LBImproved: Fast Nearest-Neighbor Retrieval under Dynamic Time Warping

[![CI](https://github.com/lemire/lbimproved/actions/workflows/ci.yml/badge.svg)](https://github.com/lemire/lbimproved/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE.txt)
[![Header-only](https://img.shields.io/badge/C%2B%2B-header--only-green.svg)](#usage)

Finding the nearest neighbor under [Dynamic Time Warping](https://en.wikipedia.org/wiki/Dynamic_time_warping) (DTW) is a core operation in time series analysis -- but DTW is expensive. A naive scan computes the full DTW distance against every candidate, which quickly becomes a bottleneck as your dataset grows.

**LBImproved** solves this with a two-pass lower-bound filter that skips most candidates *without ever computing their full DTW distance*. It is a single-header C++ library: drop one file into your project and go.

## How it works

The key insight is that a cheap lower bound on DTW distance can rule out candidates that cannot possibly beat the current best match. LB_Keogh (Keogh, 2002) provides one such bound by projecting candidates onto an envelope around the query. **LB_Improved** adds a second pass -- projecting the query onto the *candidate's* envelope -- to produce a tighter bound that prunes even more candidates.

The result: you get the **exact same nearest neighbor** as a full DTW scan, just faster.

## Quick start

The library is a single header: [`include/dtw.h`](include/dtw.h). Copy it into your project, or use CMake's `add_subdirectory`.

```cpp
#include "dtw.h"

// Your query time series
std::vector<double> query = load_query();
int constraint = query.size() / 10; // 10% warping window

// Create a filter and scan your database
LB_Improved filter(query, constraint);
for (auto& candidate : database) {
    filter.test(candidate); // internally prunes or computes full DTW
}

std::cout << "Best match distance: " << filter.getLowestCost() << std::endl;
```

That's it. The filter handles all the pruning internally. Every `test()` call either:
1. Rejects the candidate via the lower bound (fast), or
2. Computes the full DTW distance and updates the best match.

## Algorithms

All classes share the same interface: construct with a query and a warping constraint, then call `test()` on each candidate.

| Class | Strategy |
|---|---|
| `NaiveNearestNeighbor` | Full DTW on every candidate (baseline) |
| `LB_Keogh` | Single-pass lower bound pruning |
| `LB_KeoghEarly` | LB_Keogh + early termination during DTW |
| `LB_Improved` | Two-pass lower bound (tighter pruning) |
| `LB_ImprovedEarly` | LB_Improved + early termination during DTW |

Start with `LB_ImprovedEarly` for best performance. Use `NaiveNearestNeighbor` as a correctness baseline.

## Building and testing

### CMake (recommended)

```bash
cmake -B build
cmake --build build
ctest --test-dir build       # run tests
./build/benchmark            # run benchmarks
```

### Make

```bash
make
./unittesting
./benchmark
```

## Reference

Daniel Lemire, [Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound](http://arxiv.org/abs/0811.3301), *Pattern Recognition* 42(9), pages 2169-2180, 2009.

If you use this library in your research, please cite:

```bibtex
@article{Lemire2009,
  author  = {Daniel Lemire},
  title   = {Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound},
  journal = {Pattern Recognition},
  volume  = {42},
  number  = {9},
  pages   = {2169--2180},
  year    = {2009},
  doi     = {10.1016/j.patcog.2008.11.030}
}
```

Independent validation from Wang et al. (2013):

> To our knowledge, there is only one paper that offers a plausible speedup based on a tighter lower bound -- Lemire (2009) suggests a mean speedup of about 1.4 based on a tighter bound. These results are reproducible, and testing on more general data sets we obtained similar results.

## License

Apache License 2.0. See [LICENSE.txt](LICENSE.txt).

## See also

- [dtwclust](https://github.com/asardaes/dtwclust): R package for time series clustering with DTW optimizations
