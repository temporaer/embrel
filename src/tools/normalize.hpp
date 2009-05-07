#ifndef __NORMALIZE_H__
#define __NORMALIZE_H__

#include <vector>
#include "stats.hpp"

template<class I, class T=float>
class Denormalizer {
	private:
		bool mValid;
		I mBegin, mEnd;
		T mSub,mDiv,mMul,mAdd;
	public:
		Denormalizer(T sub, T div, T mul, T add, I begin, I end) :
			mValid(true),
			mBegin(begin),mEnd(end),
			mSub(sub),mDiv(div),mMul(mul),mAdd(add)
		{
		}
		void operator()(){
			(*this)(mBegin,mEnd);
		}
		void operator()(I begin, I end){
			if(!mValid) return;
			mValid = false;
			for(I i=begin; i!= end; i++)
			{
				T vi = *i - mSub;
				vi  /= mDiv;
				vi  *= mMul;
				*i   = mAdd + vi;
			}
		}
};

template<class I, class T>
inline
Denormalizer<I,T> normalize_minmax(I begin, I end, T vmin, T vmax, ExactDescriptiveStatistics* stats=NULL)
{
	T tMin,tRange;
  if(stats==NULL){
    ExactDescriptiveStatistics s(begin,end);
		tMin    = s.getMin();
		tRange  = s.getRange();
  }else{
		tMin    = stats->getMin();
		tRange  = stats->getRange();
	}
  for(I i=begin;i!=end;i++)
  {
    T vi = (*i) - tMin;
    vi  /= tRange;
    vi  *= vmax-vmin;
    (*i) = vi + vmin;
  }
	return Denormalizer<I,T>(vmin,vmax-vmin,tRange,tMin,begin,end);
}

template<class T, class V>
inline
T
normalize_sd(const T& t, const V& vmin, const V& vmax, const V& fact, const ExactDescriptiveStatistics& stats){
		T tSD = stats.getSD() * fact;
    T vi = t - stats.getMean() + tSD;
    vi  /= 2 * tSD;
    vi  *= vmax-vmin;
    vi  += vmin;
	return vi;
}
template<class T, class V>
inline
T
normalize_minmax(const T& t, const V& vmin, const V& vmax, const ExactDescriptiveStatistics& stats){
    T vi = t - stats.getMin();
    vi  /= stats.getRange();
    vi  *= vmax-vmin;
    vi  += vmin;
	return vi;
}

template<class I, class T>
inline
Denormalizer<I,T> normalize_sd(I begin, I end, T vmin, T vmax, T fact, ExactDescriptiveStatistics* stats=NULL)
{
	T tMean,tSD;
  if(stats==NULL){
    ExactDescriptiveStatistics s(begin,end);
		tMean    = s.getMean();
		tSD      = s.getSD() * fact;
  }else{
		tMean    = stats->getMean();
		tSD      = stats->getSD() * fact;
	}
	T tRange = vmax - vmin;
  for(I i=begin;i!=end;i++)
  {
		T vi = (*i) - tMean + tSD; 
    vi  /= 2 * tSD;
    vi  *= tRange;
    *i   = vi + vmin;
  }
	return Denormalizer<I,T>(vmin,tRange,2*tSD,tMean-tSD,begin,end);
}

// this "normalization" destroys information and therefore cannot be denormalized
template<class I, class T>
inline
void normalize_cutsd(I begin, I end, T fact, T vmin, T vmax, ExactDescriptiveStatistics* stats=NULL)
{
	T tMean,tSD;
  if(stats==NULL){
    ExactDescriptiveStatistics s(begin,end);
		tMean    = s.getMean();
		tSD      = s.getSD();
  }else{
		tMean    = stats->getMean();
		tSD      = stats->getSD();
	}
	T vSDRange = fact*tSD;
	T vMin     = tMean - fact*tSD;
	T vMax     = tMean - fact*tSD;
	T vRange   = vmax  - vmin;
  for(I i=begin;i!=end;i++)
  {
    T vi = max(vMin,min(*i,vMax));  // [vMin     , vMax      ]
    vi -= tMean - vSDRange;         // [0        , +2vSDRange]
    vi /= 2 * vSDRange;             // [0        , +2        ]
    vi *= vRange;                   // [0        , +2vRange  ]
    *i  = *i + vmin;                // [vmin     , vmax      ] 
  }
}

enum NormalizationType {
	NT_MINMAX,
	NT_STANDARD_DEV,
	NT_CUT_STANDARD_DEV
};


template<class I, class T>
inline
void normalize(NormalizationType nt, I begin, I end, T vmin, T vmax, T fact, ExactDescriptiveStatistics* stats=NULL)
{
	switch(nt){
		case NT_MINMAX:
			normalize_minmax(begin,end,vmin,vmax,stats);         break;
		case NT_STANDARD_DEV:
			normalize_sd(begin,end,vmin,vmax,fact,stats);        break;
		case NT_CUT_STANDARD_DEV:
			normalize_cutsd(begin,end,vmin,vmax,fact, stats);    break;
	}
}

#endif /* __NORMALIZE_H__ */
