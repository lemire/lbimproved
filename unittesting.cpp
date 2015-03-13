#include "dtw.h"

void testmaxmin() {
	vector<double> x;
	x.push_back(1);
	x.push_back(2);
	x.push_back(3);
	x.push_back(1.1);
	x.push_back(1);
	x.push_back(9);
	x.push_back(1.01);
	x.push_back(11);
	x.push_back(0.001);
	vector<double> U(x.size()),L(x.size());
	vector<double> U2(x.size()),L2(x.size());
	for(uint constraint = 0; constraint<x.size(); ++constraint) {
		
		computeEnvelope(x,constraint,U,L);
		Envelope env;
		env.compute(x,constraint,U2,L2);
		if((U!=U2) or (L!=L2)) { cout<< "bug! "<<constraint<<endl;
		cout <<" x U L U2 L2"<<endl;
		for(uint k = 0; k<x.size();++k)
		  cout<<"k="<<k<<" "<<x[k]<<" "<< U[k]<<" "<< L[k]<<" " <<U2[k]<<" "<< L2[k]<<endl;
		}
	}
	
}

void testdtw() {
	vector<double> x;
	x.push_back(1);
	x.push_back(2);
	x.push_back(3);
	x.push_back(1.1);
	x.push_back(1);
	x.push_back(9);
	x.push_back(1.01);
	x.push_back(11);
	x.push_back(0.001);
	vector<double> y;
	for(uint k = 0; k< x.size();++k)
	  y.push_back(x[k]+15);
	dtw mDTW(x.size(),1);
	double fd = mDTW.fastdynamic(x,y);
	if(fd!=l1diff(x,y)) cerr<<"bug"<<endl;
}


int main() {
	testmaxmin();
	testdtw();
	return 0;
}