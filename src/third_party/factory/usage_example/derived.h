// file: derived.h
// Copyright 2001, Jim Hyslop.
// This file may be freely used, modified and distributed, provided that
// the accompanying copyright notice remains intact.
//
// Definition of the Derived class
#ifndef DERIVED_CLASS_DEFINED
#define DERIVED_CLASS_DEFINED

#include "base.h"
#include <string>

class Derived : public Base
{
public:
    // A virtual function to illustrate how well this system
    // works.
    virtual void doSomething();
};

#endif