#include "mex.h"

#include <math.h>




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *out_prob,*prob;
	int N,i,rowin;
    int xind = 0;
       

    prob = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]);

    /*    printf("N=%d\n",N);*/
	plhs[0] = mxCreateDoubleMatrix(1,N*(N+1)/2,mxREAL);
	/*    printf("DoneAlloc");*/
	out_prob = mxGetPr(plhs[0]);

	for (i=N; i>0; i--)
	{
	  /*                printf("%d\n",i); */
		prob+=(N-i);
		memcpy(out_prob,prob,sizeof(double)*i);
		prob+=i;
		out_prob+=i;
	}
}


