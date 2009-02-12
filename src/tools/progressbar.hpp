#ifndef __PROGRESSBAR_HPP__
#define __PROGRESSBAR_HPP__
#include <iostream>
#include <string>

class ProgressBar {
  public:
    ProgressBar(int i=100, const std::string& desc = "Working" ,int cwidth=30)
      :ivpDesc(desc),
       ivMax(i),
       ivCurrent(0),
       ivClearLen(0),
       ivCWidth(cwidth),
       ivCPos(0)
    {
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
    inline void finish(){
      std::cout << "\r";
      for(uint i=0;i<ivCWidth + ivClearLen + 10 + ivpDesc.size(); i++)
	std::cout << " ";
      std::cout<<"\r"<<std::flush;
    }
    inline void finish(char* s){
      std::cout << s << std::endl<<std::flush;
    }
  private:
	std::string ivpDesc;
    int ivMax;
    int ivCurrent;
    uint ivClearLen;

    int ivCWidth;
    int ivCPos;

    void display(char* info=""){
      int newpos = (int)(ivCurrent * ivCWidth / (float) ivMax);
      if(newpos == ivCPos) return;

      ivCPos = newpos;

      std::cout << "\r"<<ivpDesc<<" |";
      for(int i=0;i<newpos;i++)
       std::cout << "#";
      for(int i=0;i<ivCWidth-newpos;i++)
       std::cout << " ";
      std::cout << "| [" << (int)(100.f*ivCurrent/ivMax) << "%] ";
      std::cout << info;
      std::cout << std::flush;
      ivClearLen = ivClearLen>strlen(info)?ivClearLen:strlen(info);
    }
};

#endif /* __PROGRESSBAR_HPP__ */
