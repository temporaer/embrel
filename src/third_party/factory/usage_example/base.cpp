// file: base.cpp
// Copyright 2001, Jim Hyslop.
// This file may be freely used, modified and distributed, provided that
// the accompanying copyright notice remains intact.

// Implementation file for class Base.

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "base.h"
#include <iostream>
#include "factory.h"

// Helper to keep maintenance down.  This is the unique
// identifier for the class.  It could just as easily be a class static
// variable.
const Base::BASE_KEY_TYPE baseName ("Base");

namespace
{
    // Note that this highlights one of the dangers of this technique.  We have
    // no control over the order in which various modules will call regCreateFn,
    // so this technique cannot be safely used if the factory is required
    // to instantiate an object before main() executes.
    registerInFactory<Base, Base> registerMe(baseName);
}

void Base::doSomething()
{
    std::cout << "Hi, I'm a Base class." << std::endl;
}

Base::~Base()
{
}

