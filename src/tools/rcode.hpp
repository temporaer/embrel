#ifndef __RCODE_HPP__
#define __RCODE_HPP__

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/shared_ptr.hpp>
#include <ublas_functors.hpp>
#include <rprop.hpp>

class RCode{
	public:
		typedef double precision;
		typedef boost::numeric::ublas::matrix<precision> mat_t;
		typedef boost::numeric::ublas::vector<precision> vec_t;
	public:
		// co-occurrence probabilities
		mat_t        mPxy;     ///< co-occurrence matrix of X and Y
		mat_t        mPxx;     ///< co-occurrence matrix of X and X

		// marginals
		vec_t        mMx;      ///< marginal of X (sum over Y in mPxy)
		vec_t        mMy;      ///< marginal of Y (sum over X in mPxy)
		vec_t        mMxx;     ///< marginal of X (sum over X in mPxx) 

		unsigned int mDim;     ///< target dimension
		mat_t        mXpos;    ///< position of X objects
		mat_t        mYpos;    ///< position of Y objects


		/// Constructor

		template <class T>
		void setPxy(const T& m, bool div=true){ 
			mPxy  = m; 
			if(div)
				mPxy /= boost::numeric::ublas::sum(mPxy); 
		}

		template <class T>
		void setPxx(const T& m, bool div=true){
			mPxx  = m; 
			if(div)
				mPxx /= boost::numeric::ublas::sum(mPxx); 
		}

		void run(int dim);
	private:
		void init_positions();
		void prepare_marginals();


		/// function pointer to the derivative of the log-likelihood 
		precision (RCode::*dl_dphix) (unsigned int x, unsigned int dim);
		/// function pointer to the derivative of the log-likelihood 
		precision (RCode::*dl_dpsiy) (unsigned int x, unsigned int dim);

		precision pcm_xgrad(unsigned int x, unsigned int dim);
		precision pcm_ygrad(unsigned int y, unsigned int dim);
};

#endif /* __RCODE_HPP__ */

