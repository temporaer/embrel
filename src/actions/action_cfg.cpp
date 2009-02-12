#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

class ActionCfg{
	public:
		ActionCfg();
};

ActionCfg::ActionCfg(){
	options_description od("========== Action Options ==========");

	options_description count("  Count Options");
	count.add_options()
		("count.recount,r", bool_switch(), "Force re-counting")
		;
	od.add(count);
	gCfg().addModuleOptions(od);
}

namespace {
	ActionCfg _actioncfg;
}
