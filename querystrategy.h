/**
* (c) Daniel Lemire, 2008
* This is a reimplementation of the R-tree based time series indexing approach,
* as developed by Keogh et al. initially.
*
* You must first build and install the spatial index library (http://research.att.com/~marioh/spatialindex/index.html)
* I built this software with release 1.3.2 - May 23rd, 2008.
*/

#ifndef QUERYSTRAG
#define QUERYSTRAG
#include <algorithm>
#include <SpatialIndex.h>

using namespace SpatialIndex;
using namespace std;


/**
* You should never need to use this class directly.
*/
class DTWQueryStrategy : public SpatialIndex::IQueryStrategy {
	public:
	DTWQueryStrategy(const vector<double> &  query, const uint constraint, const uint dim, NearestNeighbor & nn,TimeSeriesDB& d) : 
	mQuery(query), U(dim), L(dim), 
	Q(),
	 NN(nn), epsilon(NN.getLowestCost()), diskfile(d), unprunedcandidates() {
		vector<floattype> maxvalues (query.size()), minvalues(query.size());
		computeEnvelope(query, constraint, maxvalues, minvalues);
		piecewiseSumReduction(maxvalues, U);
		piecewiseSumReduction(minvalues, L);
	} 
	void getNextEntry(const IEntry& ie, id_type& nextEntryToFetch, bool& bFetchNextEntry) {
			const INode & n (reinterpret_cast<const INode &>(ie));
			if (n.getLevel() == 0) {
				for (size_t cChild = 0; cChild < n.getChildrenCount(); cChild++) {
					IShape * p; 
					n.getChildShape(cChild, &p);
					Region x;
					p->getMBR(x);
					double lb = computeLowerBound(x);// *(n.m_ptrMBR[cChild]));
					delete p;
					if(lb < epsilon) {// we have a candidate!
						if(verifynodouble) {
							assert(unprunedcandidates.find(n.getChildIdentifier(cChild))== unprunedcandidates.end());
						    unprunedcandidates.insert(n.getChildIdentifier(cChild));
						}
						NN.test(diskfile.readTimeSeries(n.getChildIdentifier(cChild)));
						epsilon = NN.getLowestCost();
					}
				}
			} else {
				for (size_t cChild = 0; cChild < n.getChildrenCount(); cChild++)	{
					// here, attempt to discard the shape if we can!
					IShape * r;
					n.getChildShape(cChild, reinterpret_cast<IShape**>(&r));
					Region x;
					r->getMBR(x);
					double lb = computeLowerBound(x);// *(n.m_ptrMBR[cChild]));
					delete r;
					if(lb<epsilon) 
						Q.push(pair<double,id_type>(lb,n.getChildIdentifier(cChild)));
				}
			}
			if(Q.empty())
				bFetchNextEntry=false;
			else {
				// here I need to loop through the stack
				 pair<double,id_type>  nextbest ( Q.top() );
				 Q.pop();
				 if(nextbest.first < epsilon) {
				 	bFetchNextEntry=true;
				 	nextEntryToFetch = nextbest.second;
				 } else {
				   bFetchNextEntry=false;
				 }
			}
	}
	
	enum{verifynodouble=false};
	
	double computeLowerBound(Region & candidate) {
		assert(candidate.getDimension()==U.size());
		double cost = 0.0;
		for (uint k = 0; k< U.size(); ++k) {
			if(candidate.getHigh(k)<L[k]) {
				cost += L[k] - candidate.getHigh(k); 
			} else if (candidate.getLow(k)>U[k]) {
				cost += candidate.getLow(k) - U[k];
			}
		}
		return cost;
	}
	
	double computeLowerBoundVerbose(Region & candidate) {
		assert(candidate.getDimension()==U.size());
		cout<<"comparing MBR:"<<endl;
		for (uint k = 0; k< U.size(); ++k) 
			cout<<candidate.getLow(k)<<" "<<candidate.getHigh(k)<<" / "<<L[k]<< " "<<U[k]<<" mQuery:"<<mQuery[k]<<endl;
		double cost = 0.0;
		for (uint k = 0; k< U.size(); ++k) {
			if(candidate.getHigh(k)<L[k]) {
				cost += L[k] - candidate.getHigh(k); 
			} else if (candidate.getLow(k)>U[k]) {
				cost += candidate.getLow(k) - U[k];
			}
		}
		cout <<"cost is "<<cost<<endl;
		return cost;
	}
	
	vector<double> mQuery, U,L;
	priority_queue<pair<double,id_type>,vector<pair<double,id_type> >, greater<pair<double,id_type>  >  > Q;
	NearestNeighbor & NN;
	double epsilon; // initialized to be very large
	TimeSeriesDB &diskfile;
	set<id_type> unprunedcandidates;

};

#endif 

