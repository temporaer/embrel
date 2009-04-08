#ifndef __UBLAS_FUNCTORS_HPP__
#define __UBLAS_FUNCTORS_HPP__

//vector<double> v( 100 );
// // ... fill v ...                                     
//std::cout << apply_to_all( v, functor::log<double>() ) << std::endl; 
//std::cout << apply_to_all<functor::log<double> >( v ) << std::endl; 
	

#include <boost/numeric/ublas/fwd.hpp>

namespace boost{
namespace numeric{
namespace ublas{

// (op v) [i] = op( v [i] )
template<class OP, class E> 
BOOST_UBLAS_INLINE
typename vector_unary_traits<E, OP>::result_type
apply_to_all (const vector_expression<E> &e, const OP& op = OP() ) {
	typedef typename vector_unary_traits<E, OP>::expression_type expression_type;
	return expression_type (e ());
}


namespace functor {
	template <class T> struct log {
				typedef T value_type;
				typedef T result_type;
				log() { }
				static result_type apply(const value_type& x) { return std::log(x); }
		};
	template <class T, class R=signed char> struct sgn {
				typedef T value_type;
				typedef R result_type;
				sgn() { }
				static result_type apply(const value_type& n) {
						if (n < 0) return -1;
						if (n > 0) return 1;
						return 0;
					}
		};
}
} // ublas
} // numeric
} // boost

#endif /* __UBLAS_FUNCTORS_HPP__ */
