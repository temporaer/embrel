#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR 1
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <fstream>
#include <numeric>
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

void RCode::configure(){
  mUse_Pxx = !gCfg().getBool("code.dont_use_pxx");
  mUse_Pxy = !gCfg().getBool("code.dont_use_pxy");
	mModel   =  gCfg().getString("code.model");
	mRPropMaxIter = gCfg().getInt("code.rprop_maxiter");
}

double
RCode::run(int globaliter){
	ofstream os((boost::format("/tmp/rcod%02d.txt")%globaliter).str().c_str());
	if(!mPosInitialized)
		init_positions();
	prepare_marginals();

	dl_dphix = &RCode::pcm_xgrad;
	dl_dpsiy = &RCode::pcm_ygrad;

	RProp rp(mPxy.size1()*mDim + mPxy.size2()*mDim);

	int iter=0;
	precision lastloglik(0);
	ProgressBar pb(mRPropMaxIter,"RCode");
	RunningDescriptiveStatistics stats(20);
	for(;iter<mRPropMaxIter;iter++){

		// calculate gradient
		precision loglik = calculate_gradient(rp);
		os   << boost::format("%03.9f")%-loglik << endl;
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
RCode::calculate_gradient(const mat_t& pxy, 
	const vec_t& mx, const vec_t& my, 
	const mat_t& xpos, const mat_t& ypos, 
	const vec_t& a, const vec_t& b, 
	const vec_t& a_mult, const vec_t& b_mult, 
	mat_t& xgrad, mat_t& ygrad)
{
	ofstream os("/tmp/rcod_a.txt");
	copy(a.begin(),a.end(),ostream_iterator<double>(os,"\n"));
	os.close();

	unsigned int nx = pxy.size1();
	unsigned int ny = pxy.size2();
	precision max_elem(-1E9);
	for(unsigned int y=0;y<ny;y++){
		matrix_row<const mat_t> psi(ypos,y);
		matrix_row<const mat_t>::iterator psi_begin=psi.begin(), psi_end=psi.end();
		for(unsigned int x=0;x<nx;x++){
			precision dist2(0);
			matrix_row<const mat_t> phi(xpos,x);
			matrix_row<const mat_t>::iterator p_phi=phi.begin(), p_psi=psi_begin;
			do{ precision d = *p_phi++-*p_psi++; dist2+=d*d; }while(p_psi!=psi_end);
			precision tmp = - dist2 + a(x)+b(y);
			if(max_elem < tmp)
				max_elem=tmp;
		}
	}
	mat_t phi_expect_m(ny,mDim,0), psi_expect_m(nx,mDim,0);
	vec_t px(nx,0), py(ny,0);
	precision loglik=0, Z=0;
	for(unsigned int y=0;y<ny;y++){
		matrix_row<const mat_t> psi(ypos,y);
		matrix_row<const mat_t>::iterator psi_begin=psi.begin(), psi_end=psi.end();
		for(unsigned int x=0;x<nx;x++){
			precision tmp, dist2(0);
			matrix_row<const mat_t> phi(xpos,x);
			matrix_row<const mat_t>::iterator p_phi=phi.begin(), p_psi=psi_begin;
			do{ precision d = *p_phi++-*p_psi++; dist2+=d*d; }while(p_psi!=psi_end);
#if 1
			tmp           = -dist2+a(x)+b(y)-max_elem; // in code_grad.c; max-elem makes sure exp(0) is largest value, regardless of current dist(x,y)
			loglik       += pxy(x,y)*tmp; // in code_grad.c
			tmp           = exp(tmp);     // in code_grad.c
#else
			// my version. does not really work as expected.
			loglik       += pxy(x,y)*(-dist2+a(x)+b(y)-max_elem); 
			tmp           = a_mult(x)*b_mult(y)*exp(-dist2);   
#endif
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
	grad *= 0; // accumulate gradient here
	precision loglik(0);

	unsigned int nx = mPxy.size1();
	unsigned int ny = mPxy.size2();
	unsigned int ystart=mDim*nx;
	mat_t xgrad, ygrad;
	vec_t a(nx),b(ny),a_mult(nx),b_mult(ny);
	switch(mModel[0]){
		case 'M': 
			noalias(a)      = apply_to_all<functor::log<double> >(const_cast<const vec_t&>(mMx)); 
			noalias(a_mult) = a;
			break;
		case 'U': 
			noalias(a) = zero_vector<double>(nx); 
			noalias(a_mult) = scalar_vector<double>(nx,1.0);
			break;
		default : throw runtime_error(string("Unknown Embedding Model")+mModel);
	}
	switch(mModel[1]){
		case 'M': 
			noalias(b)      = apply_to_all<functor::log<double> >(const_cast<const vec_t&>(mMy)); 
			noalias(b_mult) = b;
			break;
		case 'U': 
			noalias(b) = zero_vector<double>(ny); 
			noalias(b_mult) = scalar_vector<double>(ny,1.0);
			break;
		default : throw runtime_error(string("Unknown Embedding Model")+mModel);
	}

	if(mUse_Pxy){
		loglik += calculate_gradient(mPxy,mMx,mMy,mXpos,mYpos,a,b,a_mult,b_mult,xgrad,ygrad);

		if(!mFixXpos)
		for(unsigned int d=0;d<mDim;d++){
			vector_range<RProp::vec_t> xrange (grad, range(d*nx,(d+1)*nx));
			xrange += column(xgrad,d);
		}

		if(!mFixYpos)
		for(unsigned int d=0;d<mDim;d++){
			vector_range<RProp::vec_t> yrange (grad, range(ystart+d*ny,ystart+(d+1)*ny));
			yrange += column(ygrad,d);
		}
	}

	if(mUse_Pxx && !mFixXpos){
		double alpha = (double)nx/ny;
		//double alpha = (double)ny/nx;
		loglik += alpha*calculate_gradient(mPxx,mMxx,mMxx,mXpos,mXpos,a,a,a_mult,a_mult,xgrad,ygrad);
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
	random_init(mXpos,-1.0,1.0);
	random_init(mYpos,-1.0,1.0);
	mPosInitialized=true;
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
