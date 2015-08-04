/**
* (c) Daniel Lemire, 2008
* This is a reimplementation of the R-tree based time series indexing approach,
* as developed by Keogh et al. initially.
*
* You must first build and install the spatial index library (http://research.att.com/~marioh/spatialindex/index.html)
* I built this software with release 1.3.2 - May 23rd, 2008.
*/

#ifndef RTREEBASED
#define RTREEBASED
#include <algorithm>
#include <SpatialIndex.h>
#include <fstream>
#include <iostream>


class TimeSeriesDB {
	public:
	virtual vector<double> & readTimeSeries(uint myid) = 0;
};

#include "querystrategy.h"

using namespace std;


/**
* A time series database
*/
class TimeSeriesTree : public TimeSeriesDB{
	public:

	/*
	* constructor to use when the data is already present
	*/
	TimeSeriesTree( string  filename ) :	mFilename(filename),
	mConstraint(0),
	mDim(0),
     diskfile(NULL),file(NULL), tree(NULL), mNextID(0),
    mfile(NULL), timeserieslength(0),outbuffer(0),nbrdtw(0),candidates(0),nbrqueries(0) {
			if(!lock()) {
				cerr <<"could not lock database "<<mFilename<<endl;
				cerr <<"delete the lock file if you think this is an error"<<endl;
			}
    		string metadatafile (filename);
    		metadatafile.append(".meta");
    		ifstream in(metadatafile.c_str());
    		cout<<" about to open "<<metadatafile<<endl;
    		if(!in) {
    			cerr<<"You are trying to open a tree that does not exists!"<<endl;
    			throw "could not open metadata file";
    		}
    		in >> mDim;
    		in >> mConstraint;
    		in >> mNextID;
    		in.close();
			struct stat csv_stat;
			if(stat(mFilename.c_str(),&csv_stat) == 0) {
				//cout <<"opening "<<mFilename<<endl;
			} else {
				cout <<"file "<<mFilename<<" should exists!"<<endl;
			}
			mfile=std::auto_ptr<fstream>(new fstream(mFilename.c_str(),fstream::binary|fstream::in|fstream::out));
			if(!*mfile)cerr<<"problems! can't open file "<<filename<<endl;
			try {
			diskfile = std::auto_ptr<SpatialIndex::IStorageManager>(StorageManager::loadDiskStorageManager(filename));
			file = std::auto_ptr<SpatialIndex::StorageManager::IBuffer>(StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false));
			tree = std::auto_ptr<SpatialIndex::ISpatialIndex>(RTree::loadRTree(*file, 1));
			} catch (Tools::Exception& e) {
                cerr << "******ERROR******" << endl;
                std::string s = e.what();
                cerr << s << endl;
                throw "stop here";
            }
			IStatistics  * stat;
			tree->getStatistics(&stat);
			assert(mNextID == stat->getNumberOfData());
			delete stat;
	}
	
	/**
	* constructor to call the first time you create
	* the database
	*/
	TimeSeriesTree( string  filename, const uint constraint, const uint dim, const int parent_capacity = 100, const int leave_capacity = 100) : 
	mFilename(filename),
	mConstraint(constraint),
	mDim(dim),
     diskfile(NULL),file(NULL), tree(NULL), mNextID(0),
    mfile(NULL), timeserieslength(0),outbuffer(0),nbrdtw(0),candidates(0),nbrqueries(0) {
			mfile=std::auto_ptr<fstream>(new fstream(filename.c_str(), fstream::binary | fstream::trunc | fstream::out| fstream::in));
			if(!mfile->good())cerr<<"problems!"<<endl;
			diskfile = std::auto_ptr<SpatialIndex::IStorageManager>(SpatialIndex::StorageManager::createNewDiskStorageManager(filename, 4096));
			file = std::auto_ptr<SpatialIndex::StorageManager::IBuffer>(SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false));
			// Create a new, empty, RTree with dimensionality x, minimum load 70%, using "file" as
			// the StorageManager and the RSTAR splitting policy.
			id_type indexIdentifier;
			size_t indexCapacity (parent_capacity);
			size_t leafCapacity(leave_capacity);
			tree = std::auto_ptr<SpatialIndex::ISpatialIndex>(RTree::createNewRTree(*file, 0.7,  indexCapacity, leafCapacity, mDim, SpatialIndex::RTree::RV_RSTAR, indexIdentifier));
	}

	void close() {
		delete tree.release();
		delete file.release();
		delete diskfile.release();
		mfile->close();
		delete mfile.release();
		unlock();
	}
	
	bool lock() {
		struct stat csv_stat;
		string filename(mFilename);
		filename.append(".lock");
		if(stat(filename.c_str(),&csv_stat) != 0) {
			ofstream out(filename.c_str());
			if(!out) {
				cerr<<"can't open lock file "<<filename<<endl;
				return false;
			}
			out.close();
			return true;
		} else 
		  return false;
	}
	
	
	void unlock() {
    	string metadatafile (mFilename);
    	metadatafile.append(".meta");
    	ofstream out(metadatafile.c_str());
    	if(!out) cerr<<"could not open "<<metadatafile<<endl;
    	out<<mDim<<" "<<mConstraint<<" "<<mNextID<<endl;
    	out.close();
		string filename(mFilename);
		filename.append(".lock");
		::unlink(filename.c_str());
	}

	vector<double> & readTimeSeries(uint myid) {
		assert(myid<mNextID);
		outbuffer.resize(timeserieslength);
		if(!mfile->good()) {
			cerr <<" binary file not good "<<mFilename<<endl;
			mfile->close();
			mfile=std::auto_ptr<fstream>(new fstream(mFilename.c_str(), fstream::binary | fstream::app | fstream::out| fstream::in));
			if(!mfile->good()) cerr<<"still not good!"<<endl;
			//throw "this should never happen";
		}
		mfile->seekg(outbuffer.size()*sizeof(double)*myid);
		if(!mfile->good()) {
			cerr <<" binary file not good2"<<endl;
			throw "this should never happen";
		}
		mfile->read(reinterpret_cast<char*>(&outbuffer[0]),outbuffer.size()*sizeof(double));
		if(!mfile->good()) {
			cerr <<" binary file not good3"<<endl;
			throw "this should never happen";
		}
		return outbuffer;
	}
	
	~TimeSeriesTree() {
	}

	bool good() {
		return mfile->good();
	}
	
	void add(const vector<double> &  newdata) {
		if(timeserieslength == 0)
		  timeserieslength = newdata.size();
		else 
		  assert(timeserieslength == newdata.size()) ;
		vector<double> reduced(mDim);
		piecewiseSumReduction(newdata, reduced);
		Point p(&reduced[0], reduced.size());
		tree->insertData(0, 0, p,  mNextID++);
		if(!mfile->good()) {
			cerr<<"TS not ready for addition "<<  (mNextID-1) <<endl;
			throw "stop";
		}
		mfile->seekp((mNextID-1)*newdata.size()*sizeof(double));
		if(!mfile->good()) {
			cerr<<"could not position in stream to write TS "<<  (mNextID-1)*newdata.size()*sizeof(double) <<endl;
			throw "stop";
		}
		mfile->write(reinterpret_cast<const char*>(&newdata[0]),newdata.size()*sizeof(double));
		if(!mfile->good()) {
			cerr<<"could not add the ts correctly??? "<<  (mNextID-1) <<endl;
			throw "stop";
		}
	}
	
	enum{NAIVE,LB_KEOGH,LB_IMPROVED,LB_KEOGH_EARLY,LB_IMPROVED_EARLY,TREE,LINEAR};
	
	/**
	* return the l_1 distance between the query time series and the nearest neighbor in the
	* database. A slight modification of this code would return the actual nearest neighbor.
	*/
	double getNearestNeighborCost(const vector<double> &  query, uint type, uint mode = TREE) {
		++nbrqueries;
		std::auto_ptr<NearestNeighbor> nn;
		if(type==NAIVE)
		  nn = std::auto_ptr<NearestNeighbor>(new NaiveNearestNeighbor(query,mConstraint));
		else if(type == LB_KEOGH)
		  nn = std::auto_ptr<NearestNeighbor>(new LB_Keogh(query,mConstraint));
		else if(type == LB_IMPROVED)
		  nn = std::auto_ptr<NearestNeighbor>(new LB_Improved(query,mConstraint));
		else if(type == LB_KEOGH_EARLY)
		  nn = std::auto_ptr<NearestNeighbor>(new LB_KeoghEarly(query,mConstraint));
		else if(type == LB_IMPROVED_EARLY)
		  nn = std::auto_ptr<NearestNeighbor>(new LB_ImprovedEarly(query,mConstraint));
		else {
		  cerr <<"I don't know this type"<<endl; throw "stop this!";
		}
		if(mode == TREE) {
			DTWQueryStrategy qs(query, mConstraint, mDim, *nn,*this);
			tree->queryStrategy(qs);
		} else if (mode == LINEAR) {
			for(uint i = 0; i <mNextID;++i) {
				vector<double> & x =  readTimeSeries(i);
				nn->test(x);
			}
		} else {
			cerr<<"Unknown mode"<<endl;
		}
		nbrdtw += nn->getNumberOfDTW();
		candidates += nn->getNumberOfCandidates();
		double epsilon = nn->getLowestCost();
		return epsilon;
	}
	
	uint getSequentialSearchCandidates()  const  {
		return candidates;
	}
	
	uint getNumberOfDTW() const {
		return nbrdtw;
	}
	
	uint getNumberOfCandidates() const {
		return mNextID*nbrqueries;
	}
	
	void resetStatistics() {
		nbrdtw=0;
		candidates = 0;
		nbrqueries=0;
	}
	
	/**
	* return the number of time series stored.
	*/
	uint getSize() {return mNextID;}

	
	
	private:

	TimeSeriesTree( const  TimeSeriesTree & o) ;

	TimeSeriesTree & operator=(const TimeSeriesTree &o);

	string mFilename;
	uint mConstraint;
	uint mDim;
	std::auto_ptr<SpatialIndex::IStorageManager> diskfile;
	std::auto_ptr<SpatialIndex::StorageManager::IBuffer> file;
	std::auto_ptr<SpatialIndex::ISpatialIndex> tree;
	uint mNextID;
	std::auto_ptr<fstream> mfile;
	uint timeserieslength;
	vector<double> outbuffer;
	uint nbrdtw,candidates, nbrqueries;
};




#endif

