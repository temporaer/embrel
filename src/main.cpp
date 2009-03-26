#include <dlfcn.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <factory/factory.h>
#include <instantiate.hpp>
#include "configuration.hpp"
#include "action.hpp"

#include <boost/shared_ptr.hpp>
#include <nana.h>
#undef C
using namespace boost;
using namespace std;


int main(int argc, char* argv[]){
	void* handle;
	handle = dlopen("graphs/libgraphs.so",               RTLD_LAZY); 
	if (!handle) {
		cerr << "Cannot open library: " << dlerror() << '\n';
		return 1;
	}

	handle = dlopen("actions/libactions.so",             RTLD_LAZY); 
	if (!handle) {
		cerr << "Cannot open library: " << dlerror() << '\n';
		return 1;
	}

	gCfg().parsecfg(argc,argv);

	boost::shared_ptr<Action> action = instantiate<Action>("action");
	action->configure();
	(*action)();
}
