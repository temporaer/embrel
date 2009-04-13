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
, mDelta0(0.01f)
, mDeltaMax(50)
, mGrad(dim_grad,0)
, mOldGrad(dim_grad,0)
, mUpdateValue(dim_grad,mDelta0)
, mDeltaW(dim_grad,0)
, mGradSgn(dim_grad,0)
, mOldGradSgn(dim_grad,0)
, mProdSgn(dim_grad,0)
{
}

void RProp::update(ARPROP_evalres res){
	if(ARPROP_DIR_OK){
		mDeltaW *= 0.5;
		return;
	}
	update();
}
void RProp::update(){
	//cout << sum(mGrad)<<endl;
	
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
		if(s < 0){ d *= mNuMinus; /*cout<<"."<<flush;*/ continue;};
		if(s > 0){ d *= mNuPlus;  /*cout<<":"<<flush;*/ d = min(d, mDeltaMax);};
	}
	swap(mGrad,mOldGrad);
	swap(mGradSgn,mOldGradSgn);
}
