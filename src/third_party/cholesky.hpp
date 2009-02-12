/* Perform Cholesky decomposition.
 *
 * For positive definite matrices, these functions compute:
 *	- cholesky decomposition (upper or lower triangular as specified)
 *	- inverse
 *	- determinant
 % Each operation comes in an in-place form and a copying form.
 * 
 * Notes: 
 * (1) The cholesky decomposition is for symmetric matrices, therefore it does not
 * matter whether the matrix is row-major or column major. However, for row-major
 * matrices, the result format (lower/upper triangular) is reversed.
 % (2) The LAPACK function dpotrf_() stores the result in the requested lower/upper
 * triangle, but does not zero-out the opposite triangular part.
 *
 * TODOs:
 *	(1) Implement reciprocal condition estimator: dpocon_
 *	(2) Implement a solver: dpotrs_
 *
 * Tim Bailey 2005.
 */

#ifndef ULAPACK_CHOLESKY_HPP_
#define ULAPACK_CHOLESKY_HPP_

//#include "scalar_fill.hpp"
namespace ulapack {

	namespace ublas = boost::numeric::ublas;

	namespace detail {
		
		// Define a series of indexing operations for various triangular forms
		template <class T>
		class TriOp;

		template <>
		class TriOp<ublas::upper> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j <= i; }
			size_t left_idx(size_t /*i*/, size_t j) { return j; }
			size_t right_idx(size_t i, size_t /*j*/) { return i; }
		};

		template <>
		class TriOp<ublas::lower> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j <= i; }
			size_t left_idx(size_t i, size_t /*j*/) { return i; }
			size_t right_idx(size_t /*i*/, size_t j) { return j; }
		};

		template <>
		class TriOp<ublas::strict_upper> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j < i; }
			size_t left_idx(size_t /*i*/, size_t j) { return j; }
			size_t right_idx(size_t i, size_t /*j*/) { return i; }
		};

		template <>
		class TriOp<ublas::strict_lower> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j < i; }
			size_t left_idx(size_t i, size_t /*j*/) { return i; }
			size_t right_idx(size_t /*i*/, size_t j) { return j; }
		};

	} // namespace detail
    
	template<class F, class A>
	void scalar_fill(ublas::matrix<double, F, A> &m, double x)
	// Fill all m with scalar x.
	{
		double *p = &m.data()[0];
		const double * const end = p + m.size1() * m.size2();
		for (; p != end; ++p) 
			*p = x;
	}

	// Fill  triangular parts of m with scalar x.
	template<class TRI, class F, class A>
	void scalar_fill(ublas::matrix<double, F, A> &m, double x)
	{
		detail::TriOp<TRI> t;
		for (size_t i = 0; i < m.size1(); ++i)
			for (size_t j = 0; t.cmp_ij(i,j); ++j)
				m(t.left_idx(i,j), t.right_idx(i,j)) = x;
	} 

} // namespace ulapack

//#include "lapack_exception.hpp"
#include <exception>

namespace ulapack {

	class LapackError : public std::exception
	{
		int info_;
		const char *message_;
	protected:
		LapackError(const char *message, int info) : 
			 info_(info), message_(message) {}
	public:
		int get_info() const { return info_; }
		const char* what() const throw() { return message_; }
	};

	class NumericalError : public LapackError
	{
	public:
		NumericalError(const char *message, int info=0) :
		  LapackError(message, info) {}
	};

	class LogicalError : public LapackError
	{
	public:
		LogicalError(const char *message, int info=0) :
		  LapackError(message, info) {}
	};

} // namespace ulapack

//#include "errormacros.hpp"
#if defined(STRINGIZE_HELPER) || defined(STRINGIZE) || defined(ERROR_INFO)
#   error ULAPACK error macros have already been defined elsewhere 
#endif

#define STRINGIZE_HELPER(exp) #exp
#define STRINGIZE(exp) STRINGIZE_HELPER(exp)

#define ERROR_INFO(message) "ERROR: " message \
	"\nFILE: " __FILE__ "\nLINE: " STRINGIZE(__LINE__)


//#include "symmetry.hpp"
namespace ulapack {
	namespace ublas = boost::numeric::ublas;

	template<class F, class A>
	bool is_symmetric(const ublas::matrix<double, F, A> &m)
	{
		if (m.size1() != m.size2())
			return false;

		for (size_t i = 0; i < m.size1(); ++i)
			for (size_t j = i+1; j < m.size2(); ++j)
				if (m(i,j) != m(j,i))
					return false;

		return true;
	}

	template<class F, class A>
	void force_symmetry(ublas::matrix<double, F, A> &m, const bool upperToLower=true)
	// Make matrix symmetric by copying upper-to-lower (default) or lower-to-upper.
	{
		if (m.size1() != m.size2())
			throw LogicalError(ERROR_INFO("Matrix is not square"));

		if (upperToLower) {
            for (size_t i = 0; i < m.size1(); ++i)
				for (size_t j = i+1; j < m.size2(); ++j)
					m(j,i) = m(i,j);
		}
		else {
            for (size_t i = 0; i < m.size1(); ++i)
				for (size_t j = i+1; j < m.size2(); ++j)
					m(i,j) = m(j,i);
		}
	}
}

namespace ulapack {
	namespace ublas = boost::numeric::ublas;

	namespace detail {

		extern "C" 
		{
		// LAPACK function declarations

		// Cholesky factorization of a real symmetric positive definite matrix
#ifdef _MSC_VER
		void dpotrf( const char uplo, const int n, double *a, 
			const int lda, int *info );
#else
		void dpotrf_( const char *uplo, const int *n, double *a, 
			const int *lda, int *info );
#endif

		// Inverse of Cholesky factorised matrix
#ifdef _MSC_VER
		void dpotri( const char uplo, const int n, double *a, 
			const int lda, int *info );
#else
		void dpotri_( const char *uplo, const int *n, double *a, 
			const int *lda, int *info );
#endif
		// Solve a system of linear equations using the Cholesky factorization
#ifdef _MSC_VER
		void dpotrs( const char uplo, const int n, const int nrhs, 
			double *a, const int lda, double *b, const int ldb, int *info );
#else
		void dpotrs_( const char *uplo, const int *n, const int *nrhs, 
			double *a, const int *lda, double *b, const int *ldb, int *info );
#endif
		
		} // extern "C" 

		template <class A>
		char uplo_flag(const ublas::matrix<double, ublas::row_major, A> &/*m*/, const bool upper)
		// Compute "reversed" upper/lower flags for row-major matrices
		{
			return (upper) ? 'L' : 'U';
		}

		template <class A>
		char uplo_flag(const ublas::matrix<double, ublas::column_major, A> &/*m*/, const bool upper)
		// Compute normal upper/lower flags for column-major matrices
		{
			return (upper) ? 'U' : 'L';
		}

		template<class F, class A>
		int cholesky_basic_checked(ublas::matrix<double,F,A> &m, const bool upper)
		// Perform Cholesky decomposition, but do not zero out other-triangular part.
		// Returns 0 if matrix was actual positive-definite, otherwise it returns a
		// LAPACK info value.
		{
			if (m.size1() != m.size2())
				throw LogicalError(ERROR_INFO("Matrix is not square"));
			assert(is_symmetric(m)); 

			// Call LAPACK routine
			int info;
			char uplo = detail::uplo_flag(m, upper);
			int size = static_cast<int>(m.size1());
#ifdef _MSC_VER
			detail::dpotrf( uplo, size, &m.data()[0], size, &info );
#else
			detail::dpotrf_( &uplo, &size, &m.data()[0], &size, &info );
#endif

			// Check validity of result
			if (info < 0) 
				throw LogicalError(ERROR_INFO("Invalid argument"), info);

			return info;
		}

		template<class F, class A>
		void cholesky_basic(ublas::matrix<double,F,A> &m, const bool upper)
		// Perform Cholesky decomposition, but do not zero out other-triangular part.
		{
			int info = cholesky_basic_checked(m, upper);
			if (info > 0) 
				throw NumericalError(ERROR_INFO("Matrix is not positive definite"), info);
		}

		template<class F, class A>
		void zero_strict_triangular(ublas::matrix<double,F,A> &m, const bool upper)
		// Zero out strict-other triangular part of matrix.
		{
			if (upper) 
				scalar_fill<ublas::strict_lower>(m, 0.);
			else 
				scalar_fill<ublas::strict_upper>(m, 0.);
		}

	} // namespace detail

	// --------------------------------------------------------------------------
	// --------------------------- In-place algorithms --------------------------
	// --------------------------------------------------------------------------

	template<class F, class A>
	bool chol_checked_inplace(ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute Cholesky decomposition (in-place).
	// Returns false if matrix is not positive-definite.
	{
		if (detail::cholesky_basic_checked(m, upper))
			return false;
		detail::zero_strict_triangular(m, upper);
		return true;
	}

	template<class F, class A>
	void chol_inplace(ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute Cholesky decomposition (in-place).
	{
		detail::cholesky_basic(m, upper);
		detail::zero_strict_triangular(m, upper);
	}

	template<class F, class A>
	void inv_chol_inplace(ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute inverse (in-place) of a pos. def. matrix that has previously been 
	// Cholesky decomposed (ie, "m" is already a Cholesky matrix). 
	// WARNINGS: 
	//	(1) The value of "upper" must match the actual form of "m". This function does not check.
	//	(2) This function does *not* return the inverse of a Cholesky matrix.
	{
		// Call LAPACK routine
		int info;
		char uplo = detail::uplo_flag(m, upper);
		int size = static_cast<int>(m.size1());
#ifdef _MSC_VER
		detail::dpotri( uplo, size, &m.data()[0], size, &info );
#else
		detail::dpotri_( &uplo, &size, &m.data()[0], &size, &info );
#endif

		// Check validity of result
		if (info < 0) 
			throw LogicalError(ERROR_INFO("Invalid argument"), info);
		else if (info > 0) 
			throw NumericalError(ERROR_INFO("Inverse does not exist"), info);

		// Copy result to other triangular part (ie, make a symmetric inverse matrix)
		force_symmetry(m, upper);
	}

	template<class F, class A>
	void inv_pd_inplace(ublas::matrix<double,F,A> &m)
	// Compute inverse of a positive definite matrix (in-place).
	{
		detail::cholesky_basic(m, true);
		inv_chol_inplace(m, true);
	}

	template<class F, class A1, class A2>
	void solve_chol_inplace(const ublas::matrix<double,F,A1> &a, ublas::matrix<double,ublas::column_major,A2> &b, const bool upper=true)
	// Solve a system of linear equations Ax=B with the previously Cholesky factorized matrix (ie, "a" is already a Cholesky matrix) (in-place).
	{
		// Call LAPACK routine
		int info;
		char uplo = detail::uplo_flag(a, upper);
		int size = static_cast<int>(a.size1());
		int nrhs = static_cast<int>(b.size2());
#ifdef _MSC_VER
		detail::dpotrs( uplo, size, nrhs, const_cast<double*>(&a.data()[0]), size, &b.data()[0], size, &info );
#else
		detail::dpotrs_( &uplo, &size, &nrhs, const_cast<double*>(&a.data()[0]), &size, &b.data()[0], &size, &info );
#endif

		// Check validity of result
		if (info < 0) 
			throw LogicalError(ERROR_INFO("Invalid argument"), info);
	}

	template<class F, class A1, class A2>
	void solve_chol_inplace(const ublas::matrix<double,F,A1> &a, ublas::vector<double,A2> &b, const bool upper=true)
	// Solve a system of linear equations Ax=B with the previously Cholesky factorized matrix (ie, "a" is already a Cholesky matrix) (in-place).
	{
		// Call LAPACK routine
		int info;
		char uplo = detail::uplo_flag(a, upper);
		int size = static_cast<int>(a.size1());
		int nrhs = 1;
#ifdef _MSC_VER
		detail::dpotrs( uplo, size, nrhs, const_cast<double*>(&a.data()[0]), size, &b.data()[0], size, &info );
#else
		detail::dpotrs_( &uplo, &size, &nrhs, const_cast<double*>(&a.data()[0]), &size, &b.data()[0], &size, &info );
#endif

		// Check validity of result
		if (info < 0) 
			throw LogicalError(ERROR_INFO("Invalid argument"), info);
	}

	// --------------------------------------------------------------------------
	// ------------------------- Non-in-place algorithms ------------------------
	// --------------------------------------------------------------------------

	template<class F, class A>
	ublas::matrix<double,F,A> 
	chol(const ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute Cholesky decomposition, and return result.
	{
		ublas::matrix<double,F,A> c(m);
		chol_inplace(c, upper);
		return c;
	}

	template<class F, class A>
	ublas::matrix<double,F,A> inv_pd(const ublas::matrix<double,F,A> &m)
	// Compute inverse of a positive definite matrix, and return the result.
	{
		ublas::matrix<double,F,A> inv(m);
		inv_pd_inplace(inv);
		return inv;
	}

	template<class F1, class A1, class F2, class A2>
	ublas::matrix<double,F2,A2> solve_chol(const ublas::matrix<double,F1,A1> &a, const ublas::matrix<double,F2,A2> &b)
	// Solve a system of linear equations Ax=B using Cholesky factorization, and return x.
	{
		ublas::matrix<double,F1,A1> c(a);
		detail::cholesky_basic(c, true);
		ublas::matrix<double,ublas::column_major,A2> sln(b);
		solve_chol_inplace(c, sln, true);
		return sln;
	}

	template<class F, class A1, class A2>
	ublas::vector<double,A2> solve_chol(const ublas::matrix<double,F,A1> &a, const ublas::vector<double,A2> &b)
	// Solve a system of linear equations Ax=B using Cholesky factorization, and return x.
	{
		ublas::matrix<double,F,A1> c(a);
		detail::cholesky_basic(c, true);
		ublas::vector<double,A2> sln(b);
		solve_chol_inplace(c, sln, true);
		return sln;
	}

	template<class F, class A>
	double det_chol(const ublas::matrix<double,F,A> &m)
	// Compute determinant of Cholesky matrix. 
	{
		assert(m.size1() == m.size2());
		double d = 1.;
		for (size_t i=0; i < m.size1(); ++i)
			d *= m(i,i);
		return d;
	}

	template<class F, class A>
	double det_pd_sqrt(const ublas::matrix<double,F,A> &m)
	// Compute sqrt of determinant of pos. def. matrix. 
	{
		ublas::matrix<double,F,A> c(m);
		detail::cholesky_basic(c, true);
		double dsqrt = det_chol(c);
		assert(dsqrt > 0.);
		return dsqrt;
	}

	template<class F, class A>
	double det_pd(const ublas::matrix<double,F,A> &m)
	// Compute determinant of pos. def. matrix.
	{
		double d = det_pd_sqrt(m);
		return d*d;
	}

	// Wrap everything up as a class, templatised by M - the matrix type.
	template<class M>
	class Cholesky {
	public:
		Cholesky(const M &m, const bool upper=true) :
		  cholmat(m), first(true), is_upper(upper), 
			  info(detail::cholesky_basic_checked(cholmat, is_upper))
		  { }

		bool is_posdef() const { return info == 0; }

		const M& chol() const {
			check_posdef();
			if (first) { 
				first = false;
				detail::zero_strict_triangular(cholmat, is_upper);
			}
			return cholmat;
		}

		M inv() const {
			check_posdef();
			M invmat(cholmat);
			inv_chol_inplace(invmat, is_upper);
			return invmat;
		}

		double det_sqrt() const {
			check_posdef();
			return det_chol(cholmat);
		}

		template <class F, class A>
		ublas::matrix<typename M::value_type, ublas::column_major, A> solve(const ublas::matrix<typename M::value_type, F, A> &m) const {
			check_posdef();
			if (m.size1() != cholmat.size1())
				throw LogicalError(ERROR_INFO("Matrix does not have the expected number of rows"));
			ublas::matrix<typename M::value_type, ublas::column_major, A> solvemat(m);
			solve_chol_inplace(cholmat, solvemat, is_upper);
			return solvemat;
		}

		template <class A>
		ublas::vector<typename M::value_type, A> solve(const ublas::vector<typename M::value_type, A> &v) const {
			check_posdef();
			if (v.size() != cholmat.size1())
				throw LogicalError(ERROR_INFO("Vector does not have the expected number of elements"));
			ublas::vector<typename M::value_type, A> solvevec(v);
			solve_chol_inplace(cholmat, solvevec, is_upper);
			return solvevec;
		}

		double det() const {
			double d = det_sqrt();
			return d*d;
		}

	private:
		void check_posdef() const {
			if (info != 0)
				throw NumericalError(ERROR_INFO("Matrix is not positive definite"), info);
		}

		mutable M cholmat;
		mutable bool first;
		const bool is_upper;
		const int info;
	}; 

} // namespace ulapack

#endif
