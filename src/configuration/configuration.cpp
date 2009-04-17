/*       Created   :  10/03/2008 09:37:53 PM
 *       Last Change: So Feb 08 01:00  2009 CET
 */



// STL
#include <iostream>
#include <fstream>
#include <list>
using namespace std;

// Boost
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
using namespace boost::program_options;
namespace fs = boost::filesystem;

// Nana
#include <nana.h>

// header
#include <configuration.hpp>

/****************************************************************************
 * Implementation
 ****************************************************************************/
class Configuration::Configuration_Impl{
	public:
		Configuration_Impl();
		void addModuleOptions(const boost::program_options::options_description& od);
		boost::any get(const std::string& s);
		int parsecfg(int argc, char* argv[]);
		void conflicting_options(const std::string& s, const std::string& t);
		void dependent_options(const std::string& s, const std::string& t);

	private:
		void conflicting_options(const variables_map& vm, const std::string& s, const std::string& t);
		void dependent_options(const variables_map&vm, const std::string& s, const std::string& t);
		variables_map mVM;
		std::list<pair<string,string> > mConflicts;
		std::list<pair<string,string> > mDependencies;
		positional_options_description
			mPosOpt;
		options_description
			mGeneric,
			mConfig,
			mCmdLine,
			mConfigFile,
			mVisible,
			mHidden;
		bool mVerbose;
		bool mQuiet;
		void setGenericOptions();
};



Configuration::Configuration_Impl::Configuration_Impl()
	:mGeneric ("========== Generic Options =================")
	 ,mConfig ("========== General Config Options ==========")
	 ,mCmdLine("========== Commandline Options =============")
	 ,mVisible("Allowed Options")
	 ,mHidden ("Hidden Options")
	 ,mVerbose(false)
	 ,mQuiet(false)
{
	setGenericOptions();
}

void Configuration::Configuration_Impl::setGenericOptions()
{
		mPosOpt.add("action", 1);

		mGeneric.add_options()
			("help,h",   "print usage message")
			("cfg,c",     value<string>()->default_value("config.dat"), "use this as a config file")
			("action,a",  value<string>()->default_value("count"), "what to do")
			("verbose,v", bool_switch(&mVerbose),"be verbose")
			("quiet,q",   bool_switch(&mQuiet),"be quiet")
			;

		mConfig.add_options()
			//("output,o", value<string>()->default_value("out.dat"), "filename for output")
			("output-dir", value<string>()->default_value("."), "pathname for output (prepended to output)")
			//("output-format,f", value<string>()->default_value("RealXMLPrint"), "how to print output")
			;

		mCmdLine.add(mGeneric).add(mConfig).add(mHidden);
		mConfigFile.add(mConfig).add(mHidden);

		mVisible.add(mGeneric).add(mConfig);
		conflicting_options("verbose","quiet");
}

void Configuration::Configuration_Impl::addModuleOptions(const options_description& od)
{
	mConfigFile.add(od);
	mCmdLine.add(od);
	mVisible.add(od);
}

boost::any Configuration::Configuration_Impl::get(const std::string& s)
{
	return mVM[s].value();
}

int Configuration::Configuration_Impl::parsecfg(int argc, char* argv[]){
	try{
		store(
				command_line_parser(argc,argv)
				.options(mCmdLine)
				.positional(mPosOpt).run(),mVM);
        notify(mVM);

		ifstream config_file(mVM["cfg"].as<string>().c_str());
		if(config_file){
            if(mVM["verbose"].as<bool>())
                cout<<"Loading Config from ``"<<mVM["cfg"].as<string>()<<"''"<<endl;
            store(
                    parse_config_file(config_file,mConfigFile), mVM);
            notify(mVM);
		}else{
		    if(mVM["verbose"].as<bool>())
                cout<<"Config File `"<<mVM["cfg"].as<string>()<<"' could not be opened."<<endl;
		}
		if(mVM.count("help")) {
			cout << mVisible << endl;
			exit(0);
		}
		if(!mVM.count("action")) {
			throw logic_error("Action missing!");
		}

		list<pair<string,string> >::iterator it;
		for(it = mConflicts.begin(); it != mConflicts.end(); it++)
			conflicting_options(mVM, it->first,it->second);

		for(it = mDependencies.begin(); it != mDependencies.end(); it++)
			dependent_options(mVM, it->first,it->second);

	}catch(logic_error& e){
		//cerr << "Parameter Error: "<< e.what() << endl;
		throw e;
	}

	return 1;
}

void Configuration::Configuration_Impl::conflicting_options(
		const string& s, const string& t){
	I(s != t);
	mConflicts.push_back(make_pair(s,t));
}
void Configuration::Configuration_Impl::dependent_options(
		const string& s, const string& t){
	I(s != t);
	mDependencies.push_back(make_pair(s,t));
}
void Configuration::Configuration_Impl::conflicting_options(const variables_map& vm,
		const string& opt1, const string& opt2)
{
	if (vm.count(opt1) && !vm[opt1].defaulted()
			&& vm.count(opt2) && !vm[opt2].defaulted())
		throw logic_error(string("Conflicting options '")
				+ opt1 + "' and '" + opt2 + "'.");
}

void Configuration::Configuration_Impl::dependent_options(const variables_map& vm,
		const string& for_what, const string& required_option)
{
	if (vm.count(for_what) && !vm[for_what].defaulted())
		if (vm.count(required_option) == 0 || vm[required_option].defaulted())
			throw logic_error(string("Option '") + for_what
					+ "' requires option '" + required_option + "'.");
}

/****************************************************************************
 * Wrapping-Code
 ****************************************************************************/
Configuration::Configuration()
	:mImpl(new Configuration_Impl()){
}
int Configuration::parsecfg(int argc, char* argv[])
{
	return mImpl->parsecfg(argc,argv);
}
void Configuration::addModuleOptions(const options_description& od)
{
	mImpl->addModuleOptions(od);
}
void Configuration::conflicting_options(const std::string& s, const std::string& t)
{
	mImpl->conflicting_options(s,t);
}

void Configuration::dependent_options(const std::string& s, const std::string& t)
{
	mImpl->dependent_options(s,t);
}



boost::any Configuration::getAny(const std::string& s)
{
	return mImpl->get(s);
}
string Configuration::getOutputFile(const std::string& s)
{
	string file;
	try{
		file =  get<std::string>(s);
	}catch(...){
		file = s;
	}
	string path =  get<std::string>("output-dir");
	return (fs::path(path) / fs::path(file)).string();
}
string Configuration::getString(const std::string& s)
{
	return get<std::string>(s);
}
int Configuration::getInt(const std::string& s)
{
	return get<int>(s);
}
float Configuration::getFloat(const std::string& s)
{
	return get<float>(s);
}
bool Configuration::getBool(const std::string& s)
{
	return get<bool>(s);
}
//Configuration*  Configuration::mInstance = NULL;
Configuration& gCfg(){
	static Configuration cfg;
	return cfg;
}
