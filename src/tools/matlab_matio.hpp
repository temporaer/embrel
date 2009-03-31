#ifndef __MATLABMATIO_HPP__
#define __MATLABMATIO_HPP__


#include	<climits>
#include	<limits>
#include	<boost/numeric/ublas/fwd.hpp>
#include	<boost/static_assert.hpp>
#include	<iostream>
#include	<string>
#include 	<matio.h>

/*
 * Use MatIO library to write a Boost UBLAS matrix to a file
 */
template<long unsigned int DataTypeC, long unsigned int DataTypeM, class ME>
bool matlab_matrix_out_helper(std::string filen,const char* name, const ME& m){
        typedef typename ME::size_type size_type;
        const size_type size1 = m.size1 ();
        const size_type size2 = m.size2 ();

	mat_t *mat;
	matvar_t *matvar;
	int dims[2] = {size1,size2};

	mat = Mat_Open(filen.c_str(),MAT_ACC_RDWR);

	if(mat)
	{
		//Mat_VarDelete(mat, const_cast<char*>(name));
		matvar = Mat_VarCreate(const_cast<char*>(name),DataTypeC,DataTypeM,2,dims,(void*)(&(m.data()[0])),0);
		Mat_VarWrite( mat, matvar, 1);

		Mat_VarFree(matvar);
		Mat_Close(mat);
	}else
		return 1;
	return 0;
}

template <class T>
inline bool matlab_matrix_out(std::string filen,const char* name, const boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major>& m){
	//assert(false);
	BOOST_STATIC_ASSERT(sizeof(T)==0);
	return 1;
}
//template <>
//inline bool matlab_matrix_out(std::string filen,const char* name, const boost::numeric::ublas::matrix<float,boost::numeric::ublas::column_major>& m){
	//assert(std::numeric_limits<float>::digits == 16);
	//return matlab_matrix_out_helper<MAT_C_SINGLE, MAT_T_SINGLE>(filen,name,m);
//}
template <>
inline bool matlab_matrix_out(std::string filen,const char* name, const boost::numeric::ublas::matrix<int,boost::numeric::ublas::column_major>& m){
	assert(std::numeric_limits<int>::digits == 31);
	assert(std::numeric_limits<float>::is_signed );
	return matlab_matrix_out_helper<MAT_C_INT32, MAT_T_INT32>(filen,name,m);
}
template <>
inline bool matlab_matrix_out(std::string filen,const char* name, const boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major>& m){
	return matlab_matrix_out_helper<MAT_C_DOUBLE, MAT_T_DOUBLE>(filen,name,m);
}
template <>
inline bool matlab_matrix_out(std::string filen,const char* name, const boost::numeric::ublas::matrix<unsigned char,boost::numeric::ublas::column_major>& m){
	return matlab_matrix_out_helper<MAT_C_UINT8, MAT_T_UINT8>(filen,name,m);
}


#endif /* #ifndef __MATLABMATIO_HPP__ */
