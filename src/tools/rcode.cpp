#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR 1
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <configuration.hpp>
#include "stats.hpp"
#include "progressbar.hpp"
#include "rcode.hpp"

#define foreach BOOST_FOREACH 


#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>

using namespace boost::numeric::ublas;
using namespace std;

double
RCode::run(int dim){
	mDim=dim;
	init_positions();
	prepare_marginals();

	dl_dphix = &RCode::pcm_xgrad;
	dl_dpsiy = &RCode::pcm_ygrad;

	RProp rp(mPxy.size1()*mDim + mPxy.size2()*mDim);

	int max_iter=gCfg().getInt("code.rprop_maxiter"),iter=0;
	precision lastloglik(0);
	ProgressBar pb(max_iter,"RCode");
	RunningDescriptiveStatistics stats(50);
	for(;iter<max_iter;iter++){

		// calculate gradient
		precision loglik = calculate_gradient(rp);
		pb.inc((boost::format("loglik=%3.5f")%loglik).str().c_str(),1);
		stats += fabs(lastloglik-loglik);
		if(stats.getMean() < 0.000001) break;
	

		// calculate update
		rp.update((lastloglik<loglik)?RProp::ARPROP_DIR_OK:RProp::ARPROP_DIR_WRONG);
		//rp.update();
		lastloglik = loglik;

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
	cout << "Iter = " <<iter<< ", LogLik = "<<lastloglik <<endl;
	return lastloglik;
}

template<class Mat, class Vec>
void mat_mult_cols(Mat& m, const Vec& v){
	assert(m.size2()==v.size());
	for(unsigned int i=0;i<m.size2();i++){
		column(m,i) *= v(i);
	}
}
template<class Mat, class Vec>
void mat_mult_rows(Mat& m, const Vec& v){
	assert(m.size1()==v.size());
	for(unsigned int i=0;i<m.size1();i++){
		row(m,i) *= v(i);
	}
}

double 
RCode::calculate_gradient(const mat_t& pxy, const vec_t& mx, const vec_t& my, const mat_t& xpos, const mat_t& ypos, const vec_t& a, const vec_t& b, mat_t& xgrad, mat_t& ygrad){
	unsigned int nx = pxy.size1();
	unsigned int ny = pxy.size2();
	precision max_elem(-1E9);
	for(unsigned int x=0;x<nx;x++){
		matrix_row<const mat_t> phi(xpos,x);
		for(unsigned int y=0;y<ny;y++){
			precision dist(0);
			matrix_row<const mat_t> psi(ypos,y);
			matrix_row<const mat_t>::iterator p_phi=phi.begin(),phi_end=phi.end(), p_psi=psi.begin();
			do{ precision d = *p_phi++-*p_psi++; dist+=d*d; }while(p_phi!=phi_end);
			precision tmp = - dist + a(x)+b(y);
			if(max_elem < tmp)
				max_elem=tmp;
		}
	}
	mat_t phi_expect_m(ny,mDim,0), psi_expect_m(nx,mDim,0);
	vec_t px(nx,0), py(ny,0);
	precision loglik=0, Z=0;
	for(unsigned int y=0;y<ny;y++){
		matrix_row<const mat_t> psi(ypos,y);
		for(unsigned int x=0;x<nx;x++){
			precision tmp(0),dist(0);
			matrix_row<const mat_t> phi(xpos,x);
			matrix_row<const mat_t>::iterator p_phi=phi.begin(),phi_end=phi.end(), p_psi=psi.begin();
			do{ precision d = *p_phi++-*p_psi++; dist+=d*d; }while(p_phi!=phi_end);
			tmp = -dist+a(x)+b(y)-max_elem; // makes sure exp(0) is largest value, regardless of current dist(x,y)
			loglik += pxy(x,y)*tmp;
			tmp = exp(tmp);
			//tmp = a(x)*b(y)*exp(-dist);
			//tmp = exp(-dist);
			row(phi_expect_m,y) += tmp*row(xpos,x);
			row(psi_expect_m,x) += tmp*row(ypos,y);
			px(x) += tmp;
			py(y) += tmp;
			Z     += tmp;
		}
	}
	precision Z_inv = 1/Z;
	px           *= Z_inv;
	py           *= Z_inv;
	phi_expect_m *= Z_inv;
	psi_expect_m *= Z_inv;
	loglik       -= log(Z);

	xgrad = xpos; ygrad=ypos;
	mat_t phi_exp_d = prod(trans(pxy),xpos);
	mat_t psi_exp_d = prod(pxy,ypos);


	// x gradient
	px        -= mx; 
	mat_mult_rows(xgrad,px);
	psi_exp_d -= psi_expect_m;
	xgrad     += psi_exp_d;

	// y gradient
	py        -= my;
	mat_mult_rows(ygrad,py);
	phi_exp_d -= phi_expect_m;
	ygrad     += phi_exp_d;

	return loglik;
}

double 
RCode::calculate_gradient(RProp& rp){
	RProp::vec_t& grad = rp.getGrad();
	grad *= 0; // accumulate here
	precision loglik(0);

	unsigned int nx = mPxy.size1();
	unsigned int ny = mPxy.size2();
	unsigned int ystart=mDim*nx;
	mat_t xgrad, ygrad;
	string model = gCfg().getString("code.model");
	vec_t a=  model[0]=='M' ? mMx : zero_vector<double>(nx);
	vec_t b=  model[1]=='M' ? mMy : zero_vector<double>(ny);

	if(!gCfg().getBool("code.dont_use_pxy")){
		loglik += calculate_gradient(mPxy,mMx,mMy,mXpos,mYpos,a,b,xgrad,ygrad);
		for(unsigned int d=0;d<mDim;d++){
			vector_range<RProp::vec_t> xrange (grad, range(d*nx,(d+1)*nx));
			xrange += column(xgrad,d);
		}
		for(unsigned int d=0;d<mDim;d++){
			vector_range<RProp::vec_t> yrange (grad, range(ystart+d*ny,ystart+(d+1)*ny));
			yrange += column(ygrad,d);
		}
	}

	if(!gCfg().getBool("code.dont_use_pxx")){
		double alpha = (double)nx/ny;
		loglik += alpha*calculate_gradient(mPxx,mMxx,mMxx,mXpos,mXpos,a,a,xgrad,ygrad);
		xgrad  *= alpha;
		for(unsigned int d=0;d<mDim;d++){
			vector_range<RProp::vec_t> xrange (grad, range(d*nx,(d+1)*nx));
			xrange += column(xgrad,d);
		}
	}
	return loglik;
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
RCode::pcm_xgrad(unsigned int x, unsigned int dim, double Z){
	precision v(0);
#if 0
	//cout <<"."<<flush;

	precision tmp(0);
	for(unsigned int y2=0;y2<mPxy.size2();y2++){
		precision dif = ( mXpos(x,dim)-mYpos(y2,dim) );
		tmp += MARGINALS_TIMES(x,y2)*exp(-SQ_NORM(x,y2))*dif;
	}

	for(unsigned int x1=0;x1<mPxy.size1();x1++){
		for(unsigned int y1=0;y1<mPxy.size2();y1++){
			v += mPxy(x1,y1) * tmp;
		}
	}

	v/=Z;

	for(unsigned int y1=0;y1<mPxy.size2();y1++)
		v -= mPxy(x,y1)*( mXpos(x,dim)-mYpos(y1,dim) );
	
#endif
	return 2*v;
}
RCode::precision 
RCode::pcm_ygrad(unsigned int y, unsigned int dim, double Z){
	precision v(0);
	//cout <<"."<<flush;

#if 0
	precision tmp(0);
	for(unsigned int x=0;x<mPxy.size1();x++){
		precision dif = ( mXpos(x,dim)-mYpos(y,dim) );
		tmp += MARGINALS_TIMES(x,y)*exp(-SQ_NORM(x,y))*dif;
	}

	for(unsigned int x1=0;x1<mPxy.size1();x1++){
		for(unsigned int y1=0;y1<mPxy.size2();y1++){
			v += mPxy(x1,y1) * tmp;
		}
	}

	v/=Z;

	for(unsigned int x1=0;x1<mPxy.size1();x1++)
		v -= mPxy(x1,y)*( mXpos(x1,dim)-mYpos(y,dim) );
	
#endif
	return -2*v;
}
