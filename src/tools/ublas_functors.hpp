#ifndef __UBLAS_FUNCTORS_HPP__
#define __UBLAS_FUNCTORS_HPP__

//vector<double> v( 100 );
// // ... fill v ...                                     
//std::cout << apply_to_all( v, functor::log<double>() ) << std::endl; 
//std::cout << apply_to_all<functor::log<double> >( v ) << std::endl; 
	

#include <boost/numeric/ublas/fwd.hpp>
#include <cstdlib>
#include <iostream>

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

//sum of a matrix
template < typename T>
T sum(const matrix<T> & m)
{
	T s = T(0);
	for(int row=0; row<m.size1(); row++)
		for(int col=0; col<m.size2(); col++)
			s+=m(row,col);

	return s;
}

// rand matrix
template < typename T>
void random_init(matrix<T> & m, const T& vmin, const T& vmax)
{
	T fact = vmax-vmin;
	for(unsigned int row=0; row<m.size1(); row++)
		for(unsigned int col=0; col<m.size2(); col++){
			m(row,col) = drand48() * fact + vmin;
		}
}

 
//colSums
template < typename T >
vector<T> 
colSums(matrix<T> &data)
{
	unsigned int cols = data.size2();
	ublas::vector<T> sums(cols);

	for(unsigned int i=0; i<cols; i++)
		{
			 matrix_column<ublas::matrix<T> > mc (data, i);
			 T s = T(sum(mc));
			 sums(i) = s;
		}
	return sums;
}

//==============================
//rowSums
template < typename T >
vector<T> 
rowSums(matrix<T> &data)
{
	unsigned int rows = data.size1();
	vector<T> sums(rows);

	for(unsigned int i=0; i<rows; i++)
		{
			 matrix_row<matrix<T> > mc (data, i);
			 T s = T(sum(mc));
			 sums(i) = s;
		}
	return sums;
}


} // ublas
} // numeric
} // boost

#endif /* __UBLAS_FUNCTORS_HPP__ */
