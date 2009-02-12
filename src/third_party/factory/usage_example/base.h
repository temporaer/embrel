// base.h
// Copyright 2001, Jim Hyslop.
// This file may be freely used, modified and distributed, provided that
// the accompanying copyright notice remains intact.
//
// Declaration of the base class to be stored in the factory.

#ifndef BASE_HEADER_DEFINED
#define BASE_HEADER_DEFINED

#include <string>

// The class to be stored in the factory template.
class Base
{
public:
    typedef std::string BASE_KEY_TYPE;

    // A virtual function to illustrate how well this system
    // works.
    virtual void doSomething();

    virtual ~Base();
};

#endif // sentry
