#define BOOST_TEST_MODULE test_rprop
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/format.hpp>

#include <rprop.hpp>

using namespace std;
using namespace boost::numeric;
#define V(X) #X << "="<<(X)<<" "
#define D(X) #X << "="<<(boost::format("%0+2.2f")%(X))<<" "

struct Sq1{
	double f(double x){// maximize this with rprop.
		return -x*x; 
	}
	double fbar(double x){ // derivative of f
		return -2*x;
	}
};

BOOST_FIXTURE_TEST_SUITE( suiteSq1, Sq1 )

BOOST_AUTO_TEST_CASE( testSq1 ){
	RProp rp(1);
	double pos=5.0;
	for(int i=0;i<200;i++){
		rp.getGrad(0) = fbar(pos);
		//cout << V(pos)<<V(f(pos))<<V(fbar(pos))<<endl;
		rp.update();
		pos += rp.getDeltaW(0);
	}
	BOOST_CHECK_LT(fabs(pos),0.0000001);
}

BOOST_AUTO_TEST_SUITE_END()

struct Sq2{
	double f0(double x, double y){// maximize this with rprop.
		return -x*x *(y*y + 1); 
	}
	double f0bar_x(double x, double y){
		return -x * (y*y + 1);
	}
	double f0bar_y(double x, double y){
		return -x*x *(2*y);
	}
};

BOOST_FIXTURE_TEST_SUITE( suiteSq2, Sq2 )

BOOST_AUTO_TEST_CASE( testSq2 ){
	RProp rp(2);
	double x=5.0;
	double y=5.0;
	for(int i=0;i<200;i++){
		rp.getGrad(0) = f0bar_x(x,y);
		rp.getGrad(1) = f0bar_y(x,y);
		//cout << D(x)<<D(y)<<D(f0bar_x(x,y))<<D(f0bar_y(x,y))<<D(f0(x,y))<<endl;
		cout << x <<" "<< y << " " << f0(x,y)<< endl;
		rp.update();
		x += rp.getDeltaW(0);
		y += rp.getDeltaW(1);
	}
}

BOOST_AUTO_TEST_SUITE_END()
