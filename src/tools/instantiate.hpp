#include <stdexcept>
#include <string>
#include <boost/shared_ptr.hpp>
#include <factory/factory.h>
#include <configuration.hpp>

template<class T>
boost::shared_ptr<T> instantiate(const std::string& s){
	std::string obj_nam = gCfg().getString(s);
	boost::shared_ptr<T> obj ( genericFactory<T>::instance().create(obj_nam));
	if(!obj.get())
		throw std::logic_error(std::string("Supplied Object Type ``") + obj_nam + "'' does not exist");
	return obj;
}
