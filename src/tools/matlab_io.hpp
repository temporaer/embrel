#ifndef __MATLABIO_HPP__
#define __MATLABIO_HPP__


#include	<boost/numeric/ublas/fwd.hpp>
#include	<iostream>
#include	<string>

template<class E, class T, class ME>
inline void csv_matrix_out(std::basic_ostream<E, T> &os,ME& m){
        typedef typename ME::size_type size_type;
        size_type size1 = m.size1 ();
        size_type size2 = m.size2 ();
        std::basic_ostringstream<E, T, std::allocator<E> > s;
        s.flags (os.flags ());
        s.imbue (os.getloc ());
        s.precision (os.precision ());
        if (size1 > 0) {
            if (size2 > 0)
                s << "tmp,"<<m(0, 0);
            for (size_type j = 1; j < size2; ++ j)
                s << ',' << m(0, j);
        }
        for (size_type i = 1; i < size1; ++ i) {
            s << "" << std::endl ;
            if (size2 > 0)
                s <<"tmp,"<< m(i, 0);
            for (size_type j = 1; j < size2; ++ j)
                s << ',' << m(i, j);
        }
        os << s.str().c_str ();
}
template<class E, class T, class ME>
inline void matlab_matrix_out(std::basic_ostream<E, T> &os,const char* name, const ME& m){
        typedef typename ME::size_type size_type;
        size_type size1 = m.size1 ();
        size_type size2 = m.size2 ();
        std::basic_ostringstream<E, T, std::allocator<E> > s;
        s.flags (os.flags ());
        s.imbue (os.getloc ());
        s.precision (os.precision ());
        s << name << " = [ " <<std::endl;
        if (size1 > 0) {
            if (size2 > 0)
                s << m(0, 0);
            for (size_type j = 1; j < size2; ++ j)
                s << ' ' << m(0, j);
        }
        for (size_type i = 1; i < size1; ++ i) {
            s << ";" << std::endl ;
            if (size2 > 0)
                s << m(i, 0);
            for (size_type j = 1; j < size2; ++ j)
                s << ' ' << m(i, j);
        }
        s << "];" << std::endl;
        os << s.str().c_str ();
}
template<class E, class T, class VE>
// BOOST_UBLAS_INLINE This function seems to be big. So we do not let the compiler inline it.
inline void matlab_vector_out(std::basic_ostream<E, T> &os,
		const char* name,
		const VE &v) {
	typedef typename VE::size_type size_type;
	size_type size = v.size ();
	std::basic_ostringstream<E, T, std::allocator<E> > s;
	s.flags (os.flags ());
	s.imbue (os.getloc ());
	s.precision (os.precision ());
	s << name << " = [";
	if (size > 0)
		s << v (0);
	for (size_type i = 1; i < size; ++ i)
		s << ' ' << v (i);
	s << "];"<<std::endl;
	os << s.str ().c_str ();
}


#endif /* #ifndef __MATLABIO_HPP__ */
