#include <cmath>
#include "rprop.hpp"
#include "ublas_functors.hpp"

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;

#define vec_sgn(X) apply_to_all<functor::sgn<RProp::precision> >( X )

RProp::RProp(unsigned int dim_grad)
: mNuPlus(1.2f)
, mNuMinus(0.5f)
, mDelta0(0.1f)
, mDeltaMax(50)
, mDeltaMin(0)
, mGrad(dim_grad)
, mOldGrad(dim_grad)
, mUpdateValue(dim_grad)
, mDeltaW(dim_grad)
, mGradSgn(dim_grad)
, mOldGradSgn(dim_grad)
, mProdSgn(dim_grad)
{
	mGrad        = zero_vector<RProp::precision>(dim_grad);
	mOldGrad     = zero_vector<RProp::precision>(dim_grad);
	mUpdateValue = scalar_vector<RProp::precision>(dim_grad, mDelta0);
	mDeltaW      = zero_vector<RProp::precision>(dim_grad);
	mGradSgn     = scalar_vector<RProp::sign_t>(dim_grad,0);
	mOldGradSgn  = scalar_vector<RProp::sign_t>(dim_grad,0);
  mProdSgn     = scalar_vector<RProp::sign_t>(dim_grad,0);
}

void RProp::update_irprop_plus(ARPROP_evalres globalres){
	noalias(mGradSgn)       = vec_sgn(mGrad);
	noalias(mProdSgn)       = vec_sgn(element_prod(mGradSgn,mOldGradSgn)); // add - for minimization

	sign_vec_t::iterator ps = mProdSgn.begin(), pe = mProdSgn.end();
	vec_t::iterator uw      = mUpdateValue.begin();
	vec_t::iterator dw      = mDeltaW.begin();
	sign_vec_t::iterator cg = mGradSgn.begin();
	while(ps!=pe){
		sign_t&    prod_sgn = *ps++;
		precision& delta    = *uw++;
		precision& deltaw   = *dw++;
		sign_t&    grad_sgn = *cg++;
		if(prod_sgn > 0){ 
			delta  = min(delta*mNuPlus,mDeltaMax);  
			deltaw = grad_sgn * delta;
			continue;
		};
		if(prod_sgn < 0){ 
			delta  = max(delta*mNuMinus,mDeltaMin);  
			if(globalres == ARPROP_DIR_WRONG)
				deltaw = -deltaw;
			else
				deltaw = 0;
			grad_sgn = 0;
			continue;
		};
		deltaw = grad_sgn * delta;
	}
	swap(mGradSgn,mOldGradSgn);
}
void RProp::update(){
	// step 1: update delta w
	noalias(mGradSgn)       = vec_sgn(mGrad);
	noalias(mDeltaW)        = element_prod(mGradSgn,mUpdateValue); // for minimzation: add a "-"

	// step 2: update the update-value
	noalias(mProdSgn)       = vec_sgn(element_prod(mGradSgn,mOldGradSgn));
	sign_vec_t::iterator ps = mProdSgn.begin(), pe = mProdSgn.end();
	vec_t::iterator uw      = mUpdateValue.begin();
	while(ps!=pe){
		sign_t&    s = *ps++;
		precision& d = *uw++;
		if(s < 0){ d *= mNuMinus; continue;};
		if(s > 0){ d *= mNuPlus;  d = min(d, mDeltaMax);};
	}
	swap(mGrad,mOldGrad);
	swap(mGradSgn,mOldGradSgn);
}
