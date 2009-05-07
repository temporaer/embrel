#ifndef __CODE_DATA_GEN_HPP__
#define __CODE_DATA_GEN_HPP__

#include <map>
#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>

#include <kdtree++/kdtree.hpp>
#include <progressbar.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "action.hpp"


class observation {
	public:
		typedef int klass_type;
		boost::numeric::ublas::vector<double> mPos;
		std::string           mDesc;
		klass_type            mKlass;
		int                   mRunningNumber;
		template<class Archive>
		void serialize(Archive&ar, const unsigned int version){
			ar & mDesc;
			ar & mKlass;
			ar & mRunningNumber;
		}
};

struct feature {
	std::string mId;
	std::string mFromId;
	std::string mTargetId;
	boost::numeric::ublas::vector<double> mPos;
	std::map<observation::klass_type, int> mKlassCount;
	float       mComplexity;
	int         mRunningNumber;
	int         mFreq;
	double      mEntropy;
	observation::klass_type mBestKlass;
	float       mColorFact; 
	float       mSize;
	float       mSelectCrit;
	float       mFMeasure;
	float       mSelectCrit2;
	bool        mIgnore;

	template <class Archive>
	void serialize(Archive&ar, const unsigned int version){
		ar & mId;
		ar & mFromId;
		ar & mTargetId;
		ar & mKlassCount;
		ar & mComplexity;
		ar & mRunningNumber;
		ar & mFreq;
		ar & mEntropy;
		ar & mBestKlass;
		ar & mColorFact;
		ar & mSize;
		ar & mSelectCrit;
		ar & mFMeasure;
		ar & mSelectCrit;
		ar & mIgnore;
	}
};

class CoocReader : public CSVReader{
		
	public:
		typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> matrix_itype;
		typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> matrix_dtype;
		//typedef boost::numeric::ublas::matrix<unsigned char> matrix_dtype;
		typedef boost::shared_ptr<matrix_itype> matrix_pitype;
		typedef boost::shared_ptr<matrix_dtype> matrix_pdtype;
		typedef observation::klass_type        klass_type;
		struct featpos { 
			typedef double value_type; 
			feature* mFeat; 
			value_type operator[](unsigned int i)const{return mFeat->mPos[i];}
		};
		typedef KDTree::KDTree<2,featpos>       feattree_type;

	private:
		matrix_pitype mObsFeatMat;
		matrix_pdtype mFeatFeatMat;
		int          mLines;
		observation              mCurrObsDesc;
		ProgressBar              mProgress;

	public:
		feattree_type                 mFeatTree;
		std::map<klass_type, int>     mKlasses;
		std::vector<observation>      mObsDesc;
		std::vector<feature>          mFeaDesc;
	public:
		template<class Archive>
		void serialize(Archive&ar, const unsigned int version){
			ar & mKlasses;
			ar & mObsDesc;
			ar & mFeaDesc;
			ar & mObsFeatMat;
			ar & mFeatFeatMat;
		}
		CoocReader(int lines);
		CoocReader();
		virtual void read_field (const std::string& s, int lidx, int idx, int numfields);
		virtual ~CoocReader();

		void init_features();
		inline matrix_pitype getObsFeatMat(){return mObsFeatMat;}
		inline matrix_pdtype getFeatFeatMat(){return mFeatFeatMat;}

		void load_feat_pos(const std::string& fn);
		void load_obs_pos(const std::string& fn);

		template <class T> void visit_observations(T& visitor);
		template <class T> void visit_features(T& visitor);
};

template <class T> void CoocReader::visit_observations(T& visitor)
{
	for(int i=0;i<mObsFeatMat->size1();i++){
		visitor(mObsDesc[i], boost::numeric::ublas::matrix_row<matrix_itype>(*mObsFeatMat,i));
	}
}

template <class T> void CoocReader::visit_features(T& visitor)
{
	for(int i=0;i<mObsFeatMat->size2();i++){
		visitor(mFeaDesc[i], boost::numeric::ublas::matrix_column<matrix_itype>(*mObsFeatMat,i));
	}
}

class CODE_data_gen : public Action{
	public:
		virtual void configure();
		virtual ~CODE_data_gen();
		virtual void operator()(){ this->run(); }
		void run();
		void run_code(RCode& rc, CoocReader& cr, bool load_fea_pos);
		double getLogLik(RCode& rc, CoocReader& cr);
		void writeObsChildren(std::ostream& os, CoocReader&cr, int obsnr);
		void writeFeaChildren(std::ostream& os, CoocReader&cr, int feanr);
};
#endif
