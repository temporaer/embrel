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
		("code.input_file", value<string>(), "input file basename (w/o ext)")
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

		("code.linlog_graph_out",value<string>()->default_value("linlog.graph"),"write linlog graph to this file")
		("code.linlog_pos_in"   ,value<string>()->default_value("linlog.pos"),  "read linlog-generated positions from this file")
		("code.linlog_write"   ,bool_switch(), "write linlog graph")
		("code.linlog_read"    ,bool_switch(), "read linlog positions instead of using CODE ones")

		("code.aleph_queries",    value<string>()->default_value(""), "embed additional queries (from ALEPH, f.ex.)")

		("code.dvc",      bool_switch(), "write distance vs. co-occurrence statistics")
		("code.dvc_file", value<string>()->default_value("dist_vs_cooc.dat"), "file to write distance vs. co-occurrence statistics to")

		("code.entropy_emb,e",bool_switch(),"do a 2nd, entropy-based embedding after 1st")
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
	gCfg().conflicting_options("code.linlog_read", "code.linlog_write");
	gCfg().dependent_options("code.linlog_read", "code.dont_run_code");
	gCfg().dependent_options("code.linlog_write", "code.dont_run_code");
}

namespace {
	ActionCfg _actioncfg;
}
