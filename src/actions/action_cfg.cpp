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
		("code.hebb_iter", value<int>()->default_value(100), "hebb iterations")
		("code.dont_run_code,m", bool_switch(), "dont run matlab with CODE")
		("code.remove_sim,R", bool_switch(), "remove similar features")
		("code.img_width", value<float>()->default_value(750), "Image Width")
		("code.img_height",value<float>()->default_value(600), "Image Height")
		("code.delete_perc", value<float>()->default_value(0.98), "Delete this many percent of features")
		("code.size_fact", value<float>()->default_value(1.0), "Multiply Size with this")
		("code.pos_rand", value<float>()->default_value(0.0), "Random position increment uniform (pixels)")

		("code.nrestarts",value<int>()->default_value(1),"# random restarts")
		("code.dont_use_pxy",bool_switch(),"don't use observation-query statistics")
		("code.dont_use_pxx",bool_switch(),"don't use query-query statistics")
		("code.model",value<string>()->default_value("UM"),"CODE model to use")
		("code.rprop_maxiter",value<int>()->default_value(300),"Max Num of RPROP iterations")
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
