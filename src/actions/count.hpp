#ifndef __COUNT_HPP__
#define __COUNT_HPP__
#include "action.hpp"

class Count : public Action
{
  public:
    Count();
	virtual void operator()();
	virtual void configure();
    virtual ~Count();

  private:
	bool mVerbose;

};

#endif /* #ifndef __COUNT_HPP__ */
