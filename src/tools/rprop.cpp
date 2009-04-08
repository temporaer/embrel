#include <cmath>
#include "rprop.hpp"
#include "ublas_functors.hpp"

using namespace std;
using namespace boost::numeric::ublas;

#define vec_sgn(X) apply_to_all<functor::sgn<RProp::precision> >( X )

RProp::RProp(unsigned int dim_grad)
: mNuPlus(1.2f)
, mNuMinus(0.5f)
, mDelta0(0.1)
, mDeltaMax(50)
, mGrad(dim_grad,0)
, mOldGrad(dim_grad,0)
, mUpdateValue(dim_grad,mDelta0)
, mDeltaW(dim_grad,0)
, mGradSgn(dim_grad,0)
, mOldGradSgn(dim_grad,0)
{
}

void RProp::update(){
	noalias(mGradSgn)       = vec_sgn(mGrad);
	noalias(mDeltaW)        = element_prod(-mGradSgn,mUpdateValue);
	noalias(mProdSgn)       = vec_sgn(element_prod(mGradSgn,mOldGradSgn));
	sign_vec_t::iterator ps = mProdSgn.begin(), pe = mProdSgn.end();
	vec_t::iterator dw      = mDeltaW.begin();
	for(;ps!=pe;ps++,dw++){
		sign_t& s = *ps;
		if(s < 0){ *dw *= mNuMinus; continue;};
		if(s > 0){ *dw *= mNuPlus; *dw = max(*dw, mDeltaMax);};
	}
	swap(mGrad,mOldGrad);
	swap(mGradSgn,mOldGradSgn);
}
