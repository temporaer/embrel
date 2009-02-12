// main.cpp
// Copyright 2001, Jim Hyslop.
// This file may be freely used, modified and distributed, provided that
// the accompanying copyright notice remains intact.

// Example and test file to illustrate how to use the factory template.
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "base.h"
#include "factory.h"


// Code in this source file represents the code you would have in your
// implementation file, that keeps track of the objects you need to create.

void callVirtFunc(Base &b)
{
    b.doSomething();
}

int main()
{
    std::auto_ptr<Base> newBase = genericFactory<Base>::instance().create("Base");
    std::auto_ptr<Base> newDer = genericFactory<Base>::instance().create("Derived");
    callVirtFunc(*newBase);
    callVirtFunc(*newDer);
    return 0;
}

