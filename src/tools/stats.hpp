#ifndef __STAT_H__
#define __STAT_H__

// vim:ts=2:sw=2

#include <math.h>
#include <iostream>
#include <string>
#include <queue>
#include <cassert>

class RunningDescriptiveStatistics;
class ExactDescriptiveStatistics;

std::ostream& operator<<(std::ostream&, const RunningDescriptiveStatistics & );
std::ostream& operator<<(std::ostream&, const ExactDescriptiveStatistics & );

class RunningDescriptiveStatistics {
	private:
		float ivTau;
		float ivMean;
		float ivMean2;
		float ivSD;
		float ivVar;
		std::string ivDesc;
		int ivInitialized;
	public:
		RunningDescriptiveStatistics(float wl,std::string d="")
			:ivTau(1.f/wl),
			ivMean(0),
			ivMean2(0),
			ivSD(0),
			ivVar(0),
			ivDesc(d),
			ivInitialized(0)
	{ }
		void reset(){
			ivInitialized=0;
			ivMean = ivMean2 = ivSD = ivVar = 0;
		}
		inline void operator+= (float x){ this->notify(x); }

		template<class Iter>
		inline void notify(Iter begin, Iter end){
			for(Iter i=begin;i!=end;i++)
				this->notify(*i);
		}

		inline void notify(float f){
			if(f != f){
				std::cout << "WARNING: RunningDescriptiveStatistics::notify: Ignoring NaN"<<std::endl;
				return;
			}
			switch(ivInitialized){
				case 1:
					ivMean  += ivTau*(f   - ivMean);
					ivMean2 += ivTau*(f*f - ivMean2);
					ivVar    = ivMean2 - ivMean*ivMean;
					ivSD     = sqrt(ivVar);
					break;
				case 0:
					ivMean   = f;
					ivMean2  = f*f;
					ivInitialized++;
					break;
			}
		}
		inline std::string getDesc() const{ return ivDesc; }
		inline float getTau() const{ return ivTau; }
		inline float getMean()const{return ivMean;}
		inline float getMean2()const{return ivMean2;}
		inline float getSD()const{return ivSD;}
		inline float getVar()const{return ivVar;}
};

class ExactDescriptiveStatistics {
	private:
		int   ivN;
		float ivSum;
		float ivSqSum;
		float ivMin;
		float ivMax;
		std::string ivDesc;

		mutable bool  ivDataChanged;
		mutable float ivStoredMean;
		mutable float ivStoredVar;
		mutable float ivStoredSD;
	public:
		ExactDescriptiveStatistics(std::string d=""):ivDesc(d){ reset(); }
		template<class I>
		ExactDescriptiveStatistics(I begin, I end,std::string d=""):ivDesc(d){
			reset(); 
			for(I i=begin;i!=end;i++){
				this->notify(*i);
			}
		}
		inline void operator+= (float x){ this->notify(x);}
		inline int getN()const{ return ivN; }
		inline float getSum()const{ return ivSum; }
		inline float getSum2()const{ return ivSqSum; }
		inline float getMean()const{ 
			if(ivDataChanged) calcAll();
			return ivStoredMean; 
		}
		inline float getVar() const{ if(!ivDataChanged) return ivStoredVar;  
			if(ivDataChanged) calcAll();
			return ivStoredVar; 
		}
		inline float getSD()  const{ if(!ivDataChanged) return ivStoredSD;   
			if(ivDataChanged) calcAll();
			return ivStoredSD; 
		}
		inline std::string getDesc() const{ return ivDesc; }
		inline float getMin() const{ return ivMin; }
		inline float getMax() const{ return ivMax; }
		inline float getRange() const{return ivMax-ivMin;}

		template<class Iter>
		inline void notify(Iter begin, Iter end){
			for(Iter i=begin;i!=end;i++)
				this->notify(*i);
		}

		template<class Iter,class Eval>
		inline void notify(Iter begin, Iter end, Eval O){
			for(Iter i=begin;i!=end;i++)
				this->notify(O(*i));
		}

		inline void  notify(float x){
			if(x != x){
				std::cout << "WARNING: ExactDescriptiveStatistics::notify: Ignoring NaN. State: "<<(*this)<<std::endl;
				return;
			}
			ivDataChanged=true;
			ivN     ++;
			ivSum   += x;
			ivSqSum += pow(x,2);
			ivMin    = (x<ivMin)?x:ivMin;
			ivMax    = (x>ivMax)?x:ivMax;
		}
		inline void reset(){
			ivDataChanged = true;
			ivSum = ivSqSum = 0.f; ivN = 0;
			ivMax = (float)-1E12;
			ivMin = (float) 1E12;
		}
		inline float largeSampleCI(const float z, float& low, float& high) const {
			if(ivDataChanged) calcAll();
			const float n     = (float)ivN;
			const float delta = (z*sqrt(ivStoredVar / n));
			low = ivStoredMean - delta; high = ivStoredMean + delta;
			return delta;
		}
		float largeSampleCI80(float &low, float &high) const {
			return largeSampleCI(1.28, low, high); }
		float largeSampleCI90(float &low, float &high) const {
			return largeSampleCI(1.645, low, high); }
		float largeSampleCI95(float &low, float &high) const {
			return largeSampleCI(1.96, low, high);}
		float largeSampleCI99(float &low, float &high) const {
			 return largeSampleCI(2.5758, low, high);}

		void calcAll()const;

		template<class Iter>
		static bool isInInterval(Iter begin, Iter end, float left, float right, ExactDescriptiveStatistics* stats=NULL){
			float fMin,fMax;
			if(stats==NULL){
				ExactDescriptiveStatistics s(begin,end);
				fMin = s.getMin();
				fMax = s.getMax();
			}else{
				fMin = stats->getMin();
				fMax = stats->getMax();
			}
			return (left<=fMin && right >=fMax);
		}

		// can be used to test for randomness in data. Usually tried for many "lag"s and plotted and/or tested
		template<class Iter>
		static float autoCorrelation(Iter begin, Iter end, int lag, ExactDescriptiveStatistics* stats=NULL){
			float mean, var, n, run=0;
			if(stats == NULL)
			{
				ExactDescriptiveStatistics s(begin,end);
				mean = s.getMean();
				var  = s.getVar();
				n    = s.getN();
			}else{
				mean = stats->getMean();
				var  = stats->getVar();
				n    = stats->getN();
			}
			if(lag>n){
				std::cout << "ERROR: ExactDescriptiveStatistics::autoCorrelation: too small n!"<<std::endl;
				assert(false);
			};
			std::queue<float> q;
			int i=lag;
			while(--i>0) q.push(*(begin++));
			for(Iter it = begin; it!=end; it++){
				run += ((*it)-mean)*(q.front()-mean);
				q.pop();
			}
			return (run/(n-lag)) / var;
		}

		template<class Iter1, class Iter2>
		static double correlation(Iter1 beg1, Iter1 end1, Iter2 beg2, Iter2 end2, ExactDescriptiveStatistics* stats1=NULL, ExactDescriptiveStatistics* stats2= NULL){
			float sd1,sd2;
			if(stats1==NULL) sd1 = ExactDescriptiveStatistics(beg1,end1).getSD();
			else             sd1 = stats1->getSD();
			if(stats2==NULL) sd2 = ExactDescriptiveStatistics(beg2,end2).getSD();
			else             sd2 = stats2->getSD();
			return covariance(beg1,end1,beg2,end2)/(sd1*sd2);
		}

		template<class Iter1, class Iter2>
		static float covariance(Iter1 beg1, Iter1 end1, Iter2 beg2, Iter2 end2) {
			Iter1 i1 = beg1;
			Iter2 i2 = beg2;
			float sumx=*(i1++), sumy=*(i2++), sxy=0;
			int i=1;
			for(;i1!=end1,i2!=end2;i1++,i2++,i++){
				float x = *i1, y = *i2;
				sumx += x;
				sxy  += (x-sumx/(i+1))*(y-sumy/i);
				sumy += y;
			}
			assert(i1==end1);
			assert(i2==end2);
			return sxy/i;
		}

};

#endif /* __STAT_H__ */
