
#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

class GraphCfg{
	public:
		GraphCfg();
};

GraphCfg::GraphCfg(){
	options_description od("========== Graph Options ==========");
	options_description sdfr("  SDFReader Options");
	sdfr.add_options()
		("SDFReader.files",      value<string>(),                    "sdf-files, format: file:class,file:class")
		("SDFReader.fixed_size", value<int>() ->default_value(0),    "fixed size of molecules")
		("SDFReader.use_cache",  value<bool>()->default_value(true), "use cached graphs if existent")
		;
	od.add(sdfr);
	gCfg().addModuleOptions(od);
}
namespace { GraphCfg _cfg; }
