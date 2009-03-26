#include "mex.h"

#include <math.h>



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *out_prob,*prob;
	int N,i,rowin;
    int xind = 0;
       

    prob = mxGetPr(prhs[0]);
	N = *mxGetPr(prhs[1]); /*NumberOfElements(prhs[0]); */

/*	printf("N=%d\n",N);  */
	plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
	out_prob = mxGetPr(plhs[0]);

	for (i=N; i>0; i--)
	{
		out_prob+=(N-i);
		memcpy(out_prob,prob,sizeof(double)*i);
		prob+=i;
		out_prob+=i;
	}
}


