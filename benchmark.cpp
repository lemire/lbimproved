#include <cmath>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <time.h>

#include "dtw.h"



vector<double> getrandomwalk(uint size) {
    vector<double> data(size);
    data[0] = 0.0;
    for (uint k = 1; k < size; ++k)
        data[k] = (1.0 * rand() / (RAND_MAX)) - 0.5 + data[k - 1];
    return data;
}

vector<double> getcin() {
    float val;
    cin >> val;
    vector<double> v;
    while (cin) {
        v.push_back(val);
        cin >> val;
    }
    cout << "# Read " << v.size() << " data points. " << endl;
    return v;
}



template<class NN>
vector<uint> RunMe(vector<vector<double> > & collection, vector<vector<double> > & testcollection) {
    vector<uint> bestmatches;
    clock_t start, finish;
    start = clock();

    for(uint k = 0; k < testcollection.size() ; ++k ) {
        NN n(testcollection[k],testcollection[k].size() / 10);// window set at 10%, arbitrarily
        double current = n.getLowestCost();
        uint bestmatch = 0;
        for(uint z = 0; z < collection.size(); ++z) {
            double newc = n.test(collection[z]);
            if(newc < current) {// best candidate so far
              current = newc;
              bestmatch = z;
            }
        }
        bestmatches.push_back(bestmatch);
    }
    finish = clock();
    cout << "CPU time = " << static_cast<double>(finish - start) / CLOCKS_PER_SEC << endl;

    return bestmatches;
}


void runbenchmark() {
    uint N = 128;
    vector<vector<double> > collection;
    vector<vector<double> > testcollection;

    for(uint i = 0; i < 512; ++i) {
          collection.push_back(getrandomwalk(N));
          testcollection.push_back(getrandomwalk(N));
    }
    cout << "LB Keogh"<<endl;
    RunMe<LB_Keogh>(collection,testcollection);
    cout << "LB Keogh (early)"<<endl;
    RunMe<LB_KeoghEarly>(collection,testcollection);
    cout << "LB Improved"<<endl;
    RunMe<LB_Improved>(collection,testcollection);
    cout << "LB Improved (early)"<<endl;
    RunMe<LB_ImprovedEarly>(collection,testcollection);
   
}

int main() {
	runbenchmark();
	return 0;
}
