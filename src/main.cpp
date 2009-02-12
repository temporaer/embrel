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
	void* err;
	err = dlopen("graphs/libgraphs.so",               RTLD_LAZY); I(err);
	err = dlopen("actions/libactions.so",             RTLD_LAZY); I(err);

	gCfg().parsecfg(argc,argv);

	boost::shared_ptr<Action> action = instantiate<Action>("action");
	action->configure();
	(*action)();
}
