// file: derived.cpp
// Copyright 2001, Jim Hyslop.
// This file may be freely used, modified and distributed, provided that
// the accompanying copyright notice remains intact.

// Implementation file for class Derived.

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "derived.h"
#include "factory.h"
#include <iostream>

const Base::BASE_KEY_TYPE derivedName ("Derived");

namespace {
registerInFactory<Base, Derived, Base::BASE_KEY_TYPE> registerMe(derivedName);
}

// A virtual function to illustrate how well this system
// works.
void Derived::doSomething()
{
    std::cout << "Hi, I'm a Derived\n\tMy parent says:";
    Base::doSomething();
}

