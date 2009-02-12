#include "stats.hpp"
std::ostream& operator<<(std::ostream& os, ExactDescriptiveStatistics const& eds)
{
  char s[100];
  sprintf(s,"[%s%d: %2.4f %2.4f %2.4f]",eds.getDesc().c_str(),eds.getN(),eds.getMin(),eds.getMean(),eds.getMax());
  os << s ;
  return os;
}
std::ostream& operator<<(std::ostream& os, const RunningDescriptiveStatistics & rds)
{
  char s[100];
  sprintf(s,"[%s%2.4f +- %2.4f]",rds.getDesc().c_str(),rds.getMean(),rds.getSD());
  os << s ;
  return os;
}
void ExactDescriptiveStatistics::calcAll()const
{
	if(!ivN){
		ivStoredMean =0;
		ivStoredVar = ivStoredSD = 1;
	}
	else{
		ivStoredMean = ivSum/ivN;
		ivStoredVar  = ivSqSum/ivN - pow(ivSum/ivN,2);
		if(ivStoredVar>=-0.000001)
			ivStoredSD   = sqrt(fabs(ivStoredVar));
		else
			std::cout << "WARNING: ExactDescriptiveStatistics::calcAll: Variance negative! State: "<<std::endl;
		ivDataChanged = false;

		if(ivStoredMean != ivStoredMean
				||ivStoredVar  != ivStoredVar
				||ivStoredSD   != ivStoredSD){
			std::cout << "WARNING: ExactDescriptiveStatistics::calcAll: Producing NaN. State: "<<(*this)<<std::endl;
		}
	}
}

