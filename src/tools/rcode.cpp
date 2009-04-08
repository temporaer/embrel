#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include "progressbar.hpp"
#include "rcode.hpp"


#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>

using namespace boost::numeric::ublas;
using namespace std;

void
RCode::run(int dim){
	mDim=dim;
	init_positions();
	prepare_marginals();

	dl_dphix = &RCode::pcm_xgrad;
	dl_dpsiy = &RCode::pcm_ygrad;

	RProp rp(mPxy.size1()*mDim + mPxy.size2()*mDim);

	int max_iter=200;
	ProgressBar pb(max_iter,"RCode");
	for(int iter=0;iter<max_iter;iter++){
		pb.inc();
	
		// calculate gradient
		RProp::vec_t::iterator grad_pos = rp.getGrad().begin();
		for(unsigned int dim=0;dim<mDim;dim++)
			for(unsigned int x=0;x<mPxy.size1();x++)
				*grad_pos++ = (this->*dl_dphix)(x,dim);
		for(unsigned int dim=0;dim<mDim;dim++)
			for(unsigned int y=0;y<mPxy.size2();y++)
				*grad_pos++ = (this->*dl_dpsiy)(y,dim);

		// calculate update
		rp.update();

		// apply update
		unsigned int xnum = mPxy.size1();
		unsigned int ynum = mPxy.size2();
		unsigned int yval_start = mDim * xnum;
		for(unsigned int dim=0;dim<mDim;dim++){
			matrix_column<mat_t> xpos(mXpos, dim);
			matrix_column<mat_t> ypos(mYpos, dim);

			RProp::vec_t& v = rp.getDeltaW();

			const vector_range<RProp::vec_t> xrange (v, range(dim*xnum, (dim+1)*xnum));
			const vector_range<RProp::vec_t> yrange (v, range(yval_start+dim*ynum, yval_start+(dim+1)*ynum));
			xpos += xrange;
			ypos += yrange;
		}
	}
	pb.finish();
}

void 
RCode::init_positions(){
	mXpos = mat_t(mPxy.size1(),mDim);
	mYpos = mat_t(mPxy.size2(),mDim);
	random_init(mXpos,-2.0,2.0);
	random_init(mYpos,-2.0,2.0);
}

void
RCode::prepare_marginals(){
  mMx  = rowSums(mPxy);
  mMy  = colSums(mPxy);
  mMxx = rowSums(mPxx);
}

RCode::precision 
RCode::pcm_xgrad(unsigned int x, unsigned int dim){
	precision v(0);

	for(unsigned int y=0;y<mPxy.size2();y++){
		precision dif = ( mXpos(x,dim)-mYpos(y,dim) );
		v-= mPxy(x,y)    * dif;
		v+= mMx(x)*mMy(y)* dif;
	}
	
	return 0;
	return v;
}
RCode::precision 
RCode::pcm_ygrad(unsigned int y, unsigned int dim){
	precision v(0);

	for(unsigned int x=0;x<mPxy.size1();x++){
		precision dif = ( mXpos(x,dim)-mYpos(y,dim) );
		v-= mPxy(x,y)* dif;
		v+= mMx(x)   * dif;
	}
	
	return v;
}
