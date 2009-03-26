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

	options_description code ("  Code Options");
	code.add_options()
		("code.input_file", value<string>()->default_value("/tmp/erl/chains"), "input file basename (w/o ext)")
		("code.contrast_eq_neigh", value<int>()->default_value(30), "min num neigh for contrast equalization")
		("code.dont_run_code,m", bool_switch(), "dont run matlab with CODE")
		("code.img_width", value<float>()->default_value(750), "Image Width")
		("code.img_height",value<float>()->default_value(600), "Image Height")
		("code.view_sd_fact", value<float>()->default_value(0.98), "View values higher than this many StDev")
		;
	options_description count("  Count Options");
	count.add_options()
		("count.recount,r", bool_switch(), "Force re-counting")
		("count.min_freq,f", value<int>()->default_value(10), "Minimum Frequency")
		("count.max_level,l", value<int>()->default_value(30), "Maximum Depth")
		;
	od.add(code);
	od.add(count);
	gCfg().addModuleOptions(od);
}

namespace {
	ActionCfg _actioncfg;
}
