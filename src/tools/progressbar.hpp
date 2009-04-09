#ifndef __PROGRESSBAR_HPP__
#define __PROGRESSBAR_HPP__
#include "boost/date_time/posix_time/posix_time.hpp"
#include <iostream>
#include <string>

class ProgressBar {
  public:
    ProgressBar(long int i=100, const std::string& desc = "Working" ,int cwidth=30)
      :ivpDesc(desc),
       ivMax(i),
       ivCurrent(0),
       ivClearLen(0),
       ivCWidth(cwidth),
       ivCPos(-1)
    {
      ivStartTime = boost::posix_time::second_clock::local_time();
      display();
    }
    inline void inc(char*info,int v=1){
      ivCurrent += v;
      display(info);
    }
    inline void inc(int v=1){
      ivCurrent += v;
      display();
    }
    inline void finish(bool clear=false){
      if(clear){
        std::cout << "\r";
        for(uint i=0;i<std::min((unsigned int)79,ivCWidth + ivClearLen + 60 + ivpDesc.size()); i++)
	  std::cout << " ";
        std::cout<<"\r"<<std::flush;
      }
      else{
        std::cout << std::endl;
      }
    }
    inline void finish(char* s){
      std::cout << s << std::endl<<std::flush;
    }
  private:
	std::string ivpDesc;
    long int ivMax;
    long int ivCurrent;
    uint ivClearLen;
    boost::posix_time::ptime ivStartTime;

    int ivCWidth;
    int ivCPos;

    void display(char* info=""){
      double newpos_f = ((double)ivCurrent / (double) ivMax);
      int    newpos   = (int)(ivCWidth * newpos_f);
      if(newpos == ivCPos) return;

      ivCPos = newpos;

      std::cout << "\r"<<ivpDesc<<" |";
      for(int i=0;i<newpos;i++)
       std::cout << "#";
      for(int i=0;i<ivCWidth-newpos;i++)
       std::cout << " ";
      std::cout << "| [" << (int)(100.0*newpos_f) << "% (";
      boost::posix_time::time_duration td=boost::posix_time::time_period(ivStartTime, boost::posix_time::second_clock::local_time()).length();
      int newsec = (1.0f-newpos_f)/newpos_f * td.total_seconds();
      td = boost::posix_time::seconds(newsec);
      boost::posix_time::ptime done(boost::posix_time::second_clock::local_time()+td);
      std::cout << done <<")] ";
      std::cout << info;
      std::cout << std::flush;
      ivClearLen = ivClearLen>strlen(info)?ivClearLen:strlen(info);
    }
};

#endif /* __PROGRESSBAR_HPP__ */
