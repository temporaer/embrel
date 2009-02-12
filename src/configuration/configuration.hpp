/*       Created   :  10/03/2008 09:40:20 PM
 *       Last Change: Sat Nov 22 04:00 PM 2008 CET
 */
#ifndef __CONFIGURATION_HPP__
#define __CONFIGURATION_HPP__

#include <string>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

// forward-declare any, options_description
namespace std{
	class invalid_argument;
}
namespace boost{
	namespace program_options{
		class options_description;
	}
}

//! Configures program via config file and command line
class Configuration
{
	// the public interface is doubled by Configuration_Impl,
	// where the real work is done.
	private:
		struct Configuration_Impl;                    //!< the actual implementation is in here
		boost::shared_ptr<Configuration_Impl> mImpl;  //!< use this to access implementation
	public:
		Configuration();

		/*! \brief if you write a module, send the module-options to Configuration via this method
		 * \param od the options for your module */
		void addModuleOptions(const boost::program_options::options_description& od);

		//! \brief parse commandline parameters
		int parsecfg(int argc, char* argv[]);

		/*! \brief tell Configuration about conflicting options
		 * \param s conflicting option
		 * \param t conflicting option */
		void conflicting_options(const std::string& s, const std::string& t);

		/*! \brief tell Configuration about dependent options
		 * \param s neccessary for t
		 * \param t valid only with s */
		void dependent_options(const std::string& s, const std::string& t);

		/*! \brief get any type of parameter, extract it yourself */
		boost::any getAny(const std::string& s);

		/*! \brief get any type of parameter */
		template<class T>
		T get(const std::string& s){
			try{
				T t = boost::any_cast<T>(getAny(s));
				return t;
			}
			catch(const boost::bad_any_cast &) {
				throw std::invalid_argument(
						std::string("Parameter ``")
						.append(s)
						.append("'' not given or incorrectly referenced"));
			}
		}

		//! convenience function: get output file (prepended with path)
		std::string getOutputFile(const std::string& s);
		//! convenience function: get parameter as string
		std::string getString(const std::string& s);
		//! convenience function: get parameter as float
		float getFloat(const std::string& s);
		//! convenience function: get parameter as int
		int   getInt(const std::string& s);
		//! convenience function: get parameter as bool
		bool  getBool(const std::string& s);
};

/*! the global Configuration object. 
 * We wrapped it inside a function so it is definitely constructed before use. */
extern Configuration&  gCfg(); 

#endif /* __CONFIGURATION_HPP__ */
