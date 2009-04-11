#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR 1
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include "stats.hpp"
#include "progressbar.hpp"
#include "rcode.hpp"

#define foreach BOOST_FOREACH 

#define SQ_NORM(X,Y) \
  (pow(mXpos(X,0)-mYpos(Y,0),2.0) + \
  pow(mXpos(X,1)-mYpos(Y,1),2.0))

#define P_UM
#ifdef P_UU
#	define MARGINALS_PLUS(X,Y) 0
#	define MARGINALS_TIMES(X,Y) 1
#elif defined P_UM
#	define MARGINALS_PLUS(X,Y)  mMy(Y)
#	define MARGINALS_TIMES(X,Y) mMy(Y)
#elif defined P_MU
#	define MARGINALS_PLUS(X,Y)  mMx(X)
#	define MARGINALS_TIMES(X,Y) mMx(X)
#elif defined P_MM
#	define MARGINALS_PLUS(X,Y)  mMx(X) + mMy(Y)
#	define MARGINALS_TIMES(X,Y) mMx(X) * mMy(Y)
#endif

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

	int max_iter=300,iter=0;
	precision lastloglik(0);
	ProgressBar pb(max_iter,"RCode");
	RunningDescriptiveStatistics stats(20);
	for(;iter<max_iter;iter++){

		// calculate gradient
		precision loglik = calculate_gradient(rp);
		pb.inc((boost::format("loglik=%3.5f")%loglik).str().c_str(),1);
		stats += fabs(lastloglik-loglik);
		if(stats.getMean() < 0.00001) break;
	

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
RCode::calculate_gradient(RProp& rp){
#if 0
	precision Z(0),loglik(0);
	for(unsigned int x=0;x<mPxy.size1();x++){
		for(unsigned int y=0;y<mPxy.size2();y++){
			precision sq_sum(0);
			for(unsigned int d=0;d<mDim;d++){
				sq_sum += pow(mXpos(x,d)-mYpos(y,d),2);
			}
			Z+=MARGINALS_TIMES(x,y)*exp(-sq_sum);
		}
	}
	precision logZ(log(Z));
	ExactDescriptiveStatistics elemstat;
	for(unsigned int x=0;x<mPxy.size1();x++){
		for(unsigned int y=0;y<mPxy.size2();y++){
			precision d2 = SQ_NORM(x,y);
			loglik   += mPxy(x,y)*(-logZ-d2 + log(MARGINALS_TIMES(x,y)));
			elemstat += -d2+MARGINALS_PLUS(x,y);
		}
	}
	// my method. works, but slow.
	RProp::vec_t::iterator grad_pos = rp.getGrad().begin();
	for(unsigned int dim=0;dim<mDim;dim++)
		for(unsigned int x=0;x<mPxy.size1();x++)
			*grad_pos++ = (this->*dl_dphix)(x,dim,Z);
	for(unsigned int dim=0;dim<mDim;dim++)
		for(unsigned int y=0;y<mPxy.size2();y++)
			*grad_pos++ = (this->*dl_dpsiy)(y,dim,Z);
	return loglik;
#else
	precision max_elem(-1E9);
	for(unsigned int x=0;x<mPxy.size1();x++){
		for(unsigned int y=0;y<mPxy.size2();y++){
			precision tmp = - SQ_NORM(x,y) + MARGINALS_PLUS(x,y);
			if(max_elem < tmp)
				max_elem=tmp;
		}
	}
	unsigned int nx = mPxy.size1();
	unsigned int ny = mPxy.size2();
	mat_t phi_expect_m(ny,mDim,0), psi_expect_m(nx,mDim,0);
	vec_t px(nx,0), py(ny,0);
	precision loglik=0, Z=0;
	for(unsigned int y=0;y<ny;y++){
		matrix_row<mat_t> psi(mYpos,y);
		for(unsigned int x=0;x<nx;x++){
			precision tmp(0),dist(0);
			matrix_row<mat_t> phi(mXpos,x);
			matrix_row<mat_t>::iterator p_phi=phi.begin(),phi_end=phi.end(), p_psi=psi.begin();
			do{ precision d = *p_phi++-*p_psi++; dist+=d*d; }while(p_phi!=phi_end);
			tmp = -dist+MARGINALS_PLUS(x,y)-max_elem; // makes sure exp(0) is largest value, regardless of current dist(x,y)
			loglik += mPxy(x,y)*tmp;
			//tmp = exp(tmp);
			tmp = MARGINALS_TIMES(x,y)*exp(-dist);
			row(phi_expect_m,y) += tmp*row(mXpos,x);
			row(psi_expect_m,x) += tmp*row(mYpos,y);
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

	mat_t xgrad(mXpos), ygrad(mYpos);
	mat_t phi_exp_d = prod(trans(mPxy),mXpos);
	mat_t psi_exp_d = prod(mPxy,mYpos);
	RProp::vec_t& grad = rp.getGrad();
	unsigned int ystart=mDim*nx;


	// x gradient
	px        -= mMx; 
	mat_mult_rows(xgrad,px);
	psi_exp_d -= psi_expect_m;
	xgrad     += psi_exp_d;
	for(unsigned int d=0;d<mDim;d++){
		vector_range<RProp::vec_t> xrange (grad, range(d*nx,(d+1)*nx));
		xrange = column(xgrad,d);
	}

	// y gradient
	py        -= mMy;
	mat_mult_rows(ygrad,py);
	phi_exp_d -= phi_expect_m;
	ygrad     += phi_exp_d;
	for(unsigned int d=0;d<mDim;d++){
		vector_range<RProp::vec_t> yrange (grad, range(ystart+d*ny,ystart+(d+1)*ny));
		yrange = column(ygrad,d);
	}

	return loglik;
#endif
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
	
	return 2*v;
}
RCode::precision 
RCode::pcm_ygrad(unsigned int y, unsigned int dim, double Z){
	precision v(0);
	//cout <<"."<<flush;

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
	
	return -2*v;
}
