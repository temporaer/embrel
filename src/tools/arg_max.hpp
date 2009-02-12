/*       Created   :  18/09/08 11:08:33
 *       Last Change: Thu Sep 18 08:00 PM 2008 CEST
 */


#include <boost/function.hpp>
#include <limits>

namespace util{
	using namespace ::std;
	template<class T, class Iter, class Func> 
		Iter
		arg_min(Iter begin, Iter end, T& arg, Func func){
			Iter best  = end;
			arg        = numeric_limits<T>::max();
			while(begin!=end){
				T d = func(*begin);
				if(d < arg){
					best  = begin;
					arg   = d;
				}
				begin++;
			}
			return best;
		}
	template<class T, class Iter, class Func> 
		Iter
		arg_max(Iter begin, Iter end, T& arg, Func func){
			Iter best  = end;
			arg        = numeric_limits<T>::min();
			while(begin!=end){
				T d = func(*begin);
				if(d > arg ){
					best  = begin;
					arg   = d;
				}
				begin++;
			}
			return best;
		}

};
