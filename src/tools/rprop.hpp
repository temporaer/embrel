#ifndef __RPROP_HPP__
#define __RPROP_HPP__

#include <boost/numeric/ublas/vector.hpp>

/// RProp class takes care of gradients
class RProp{
	public:
		typedef double precision;
		typedef signed char sign_t;
		typedef boost::numeric::ublas::vector<precision> vec_t;
		typedef boost::numeric::ublas::vector<sign_t> sign_vec_t;
	private:
		// parameters
		precision   mNuPlus;    // 0 < NuMinus < 1 < NuPlus
		precision   mNuMinus;   // 0 < NuMinus < 1 < NuPlus
		precision   mDelta0;    ///< initial delta
		precision   mDeltaMax;  ///< maximum delta

		// data
		vec_t   mGrad;        ///< current gradient
		vec_t   mOldGrad;     ///< last gradient
		vec_t   mUpdateValue; ///< update values (change of mDeltaW)
		vec_t   mDeltaW;      ///< step width
		sign_vec_t mGradSgn;     ///< a vector for keeping signs
		sign_vec_t mOldGradSgn;  ///< a vector for keeping signs
		sign_vec_t mProdSgn;   ///< a vector for keeping signs
		

	public:
		RProp(unsigned int dim_grad);

		/// get and set gradient
		inline vec_t&       getGrad()     { return mGrad; }
		/// get gradient
		inline const vec_t& getGrad()const{ return mGrad; }

		/// get and set gradient
		inline precision&         getGrad(unsigned int i)     { return mGrad(i); }
		/// get gradient
		inline const precision&   getGrad(unsigned int i)const{ return mGrad(i); }

		/// get update value
		inline const vec_t& getUpdateValue()const{ return mUpdateValue; }
		/// get update value
		inline const precision&   getUpdateValue(unsigned int i)const{ return mUpdateValue(i); }

		/// get weight update
		inline const vec_t& getDeltaW()const{ return mDeltaW; }
		/// get weight update
		inline const precision&   getDeltaW(unsigned int i)const{ return mDeltaW(i); }

		/// perform a 2-step update: 1) delta-w is updated 2) update-values are updated
		void update();
};


#endif /* __RPROP_HPP__ */
