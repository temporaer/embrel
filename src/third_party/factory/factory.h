// factory.h
// Generic abstract factory template.
// Copyright 2001, Jim Hyslop.
// This file may be freely used, modified and distributed, provided that
// the accompanying copyright notice remains intact.
//
// The generic abstract factory template is an implementation of the
// Abstract Class Factory pattern, as a template (see "Design Patterns:
// Elements of Reusable Object-Oriented Software", E. Gamma, R. Helm,
// R. Johnson, J. Vlissides, Addison Wesley [1995] )
//
// To use the template, you need to provide a base class and (optionally)
// a key class. The base class must provide a unique identifier
//
// The key class must be able to be used as a key in a std::map, i.e. it must
//   - implement copy and assignment semantics
//   - provide bool operator< () const;
// Default is std::string.
//
// Steps to using the factory:
// 1 - Create the base class and its derivatives
// 2 - Register each class in the factory by instantiating a
//     registerInFactory<> template class - do this in one file only (the
//     class implementation file is the perfect place for this)
// 3 - create the object by calling create() and passing it the same
//     value used when you instantiated the registerInFactory object.
// For example:
//   base header:
//   class Base { /* whatever (don't forget the virtual dtor! */ };
//
//   base implementation:
//   registerInFactory<Base, Base, std::string> registerBase("Base");
//
//   derived header:
//   class Derived : public Base { /* whatever */ };
//
//   derived implementation:
//   registerInFactory<Base, Derived, std::string> registerDer("Derived");
//
//   code that instantiates the classes:
//   std::auto_ptr<Base> newBase = genericFactory<Base>::instance().create("Base");
//   std::auto_ptr<Base> newDerived = genericFactory<Base>::instance().create("Derived");
//
// New derivatives can be added without affecting the existing code.

#ifndef FACTORY_HEADER_DEFINED
#define FACTORY_HEADER_DEFINED
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <iostream>

typedef std::string defaultIDKeyType;


// The abstract factory itself.
// Implemented using the Singleton pattern

template <class manufacturedObj, typename classIDKey=defaultIDKeyType>
class genericFactory
{
    // a BASE_CREATE_FN is a function that takes no parameters
    // and returns an auto_ptr to a manufactuedObj.  Note that
    // we use no parameters, but you could add them
    // easily enough to allow overloaded ctors, e.g.:
    //   typedef std::auto_ptr<manufacturedObj> (*BASE_CREATE_FN)(int);
    typedef std::auto_ptr<manufacturedObj> (*BASE_CREATE_FN)();

    // FN_REGISTRY is the registry of all the BASE_CREATE_FN
    // pointers registered.  Functions are registered using the
    // regCreateFn member function (see below).
    typedef std::map<classIDKey, BASE_CREATE_FN> FN_REGISTRY;
    FN_REGISTRY registry;

    // Singleton implementation - private ctor & copying, with
    // no implementation on the copying.
	genericFactory();
	genericFactory(const genericFactory&); // Not implemented
	genericFactory &operator=(const genericFactory&); // Not implemented
public:

    // Singleton access.
    static genericFactory &instance();

    // Classes derived from manufacturedObj call this function once
    // per program to register the class ID key, and a pointer to
    // the function that creates the class.
    void regCreateFn(const classIDKey &, BASE_CREATE_FN);

    // Create a new class of the type specified by className.
    std::auto_ptr<manufacturedObj> create(const classIDKey &className) const;

};

////////////////////////////////////////////////////////////////////////
// Implementation details.  If no comments appear, then I presume
// the implementation is self-explanatory.

template <class manufacturedObj, typename classIDKey>
genericFactory<manufacturedObj, classIDKey>::genericFactory()
{
}

template <class manufacturedObj, typename classIDKey>
genericFactory<manufacturedObj, classIDKey> &genericFactory<manufacturedObj, classIDKey>::instance()
{
    // Note that this is not thread-safe!
    static genericFactory theInstance;
    return theInstance;
}

// Register the creation function.  This simply associates the classIDKey
// with the function used to create the class.  The return value is a dummy
// value, which is used to allow static initialization of the registry.
// See example implementations in base.cpp and derived.cpp
template <class manufacturedObj, typename classIDKey>
void genericFactory<manufacturedObj, classIDKey>::regCreateFn(const classIDKey &clName, BASE_CREATE_FN func)
{
    registry[clName]=func;
}

// The create function simple looks up the class ID, and if it's in the list,
// the statement "(*i).second();" calls the function.
template <class manufacturedObj, typename classIDKey>
std::auto_ptr<manufacturedObj> genericFactory<manufacturedObj, classIDKey>::create(const classIDKey &className) const
{
    std::auto_ptr<manufacturedObj> ret(0);

	// TODO: remove
	//std::cout << "creating a "<<className<<std::endl;
	//typename genericFactory<manufacturedObj,classIDKey>::FN_REGISTRY::const_iterator it;
	//for(it=registry.begin();it!=registry.end();it++)
		//std::cout << it->first<<std::endl;

    typename genericFactory<manufacturedObj,classIDKey>::FN_REGISTRY::const_iterator regEntry=registry.find(className);
    if (regEntry != registry.end()) {
        ret=(*regEntry).second();
    }
    return ret;
}

// Helper template to make registration painless and simple.
template <class ancestorType,
          class manufacturedObj,
          typename classIDKey=defaultIDKeyType>
class registerInFactory
{
public:
    static std::auto_ptr<ancestorType> createInstance()
    {
        return std::auto_ptr<ancestorType>(new manufacturedObj);
    }
    registerInFactory(const classIDKey &id)
    {
        genericFactory<ancestorType>::instance().regCreateFn(id, createInstance);
		//std::cout << "registring: " << id << std::endl;
    }
};
#endif
