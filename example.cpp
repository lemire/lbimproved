

#include "dtw.h"

// this just generates some random array
vector<double> getrandomwalk(uint size) {
  vector<double> data(size);
  data[0] = 0.0;
  for (uint k = 1; k < size; ++k)
    data[k] = (1.0 * rand() / (RAND_MAX)) - 0.5 + data[k - 1];
  return data;
}

void demo(uint size) {
  std::cout << " I generated a random walk and I will try to match it with "
               "other random walks. "
            << std::endl;

  vector<double> target = getrandomwalk(size); // this is our target
  LB_Improved filter(
      target, size / 10); // we use the DTW with a tolerance of 10% (size/10)
  double bestsofar = filter.getLowestCost();
  uint howmany = 5000;
  for (uint i = 0; i < 5000; i++) {
    vector<double> candidate = getrandomwalk(size);
    double newbest = filter.test(candidate);
    if (newbest < bestsofar) {
      std::cout << " we found a new nearest neighbor, distance (L1 norm) = "
                << newbest << std::endl;
      bestsofar = newbest;
    }
  }
  std::cout << " I compared it with " << howmany
            << " random walks, closest match is at a distance (L1 norm) of "
            << filter.getLowestCost() << std::endl;
}

int main() { demo(10000); }
