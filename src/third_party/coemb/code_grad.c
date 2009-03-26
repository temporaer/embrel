#include "mex.h"
#include "math.h"

#define GetMatInd(a,xi,yi) (yi*a.nx+xi)
#define GetMatEl(a,xi,yi)  (a.dat[GetMatInd(a,xi,yi)])

double max_element;

typedef struct mat
{
	double *dat;
	int nx,ny,nall;
} Matrix;

int bDbg=1;

double dot(double *a, double *b, int n)
{
  double s=0;
  double *ap,*bp;
  double *aend=&a[n];
  for (ap=&a[0],bp=&b[0]; ap<aend; ap++, bp++)
      s+=*ap*(*bp);
  return s;
}

double vec_sum(double *a, int n)
{
  double s=0;
  double *ap;
  double *aend=&a[n];
  for (ap=&a[0]; ap<aend; ap++)
      s+=*ap;
  return s;
}

void vec_mul_scalar(double *a, int n, double sc)
{
  double s=0;
  double *ap;
  double *aend=&a[n];
  for (ap=&a[0]; ap<aend; ap++)
      (*ap)*=sc;
}

void vec_subtract(double *a, double *b, int n)
{
  double s=0;
  double *ap;
  double *bp = &b[0];
  double *aend=&a[n];
  for (ap=&a[0]; ap<aend; ap++,bp++)
      (*ap)=(*ap-*bp);
}

double l2_dist(double *a, double *b, int n)
{
  double s=0;
  double *ap,*bp;
  double *aend=&a[n];
  double d=0;

  for (ap=&a[0],bp=&b[0]; ap<aend; ap++, bp++)
  {
	  d = *ap-*bp;
      s+=(d*d);
  }
  return s;
}

void mat_mul(Matrix A, Matrix B, Matrix C)
{
	int xi,yi,j;
	double tmp;
	double *c_ptr = C.dat;
	double *b_col_ptr; 	/* Point to the beginning of the yi'th column in B */
	double *tmp_b_ptr;  /* For going over the column */
	
	b_col_ptr= B.dat;
	for (yi=0; yi<B.ny; yi++)
	{
		for (xi=0; xi<A.nx; xi++)
		{
			tmp=0;
			tmp_b_ptr = b_col_ptr;		
			for (j=0; j<A.ny; j++,tmp_b_ptr++)
				tmp+= GetMatEl(A,xi,j)*(*tmp_b_ptr);
/*				tmp+= GetMatEl(A,xi,j)*GetMatEl(B,j,yi); */
			(*c_ptr) = tmp;
			c_ptr++;
/*			C.dat[GetMatInd(C,xi,yi)] = tmp; */
		}
		b_col_ptr+= B.nx;
	}
}

/* Do res = m1'*m2; */
void mat_mul_a_tr_b(Matrix m1, Matrix m2, Matrix res)
{
	int i1,i2,k;
	double *p1,*p2;
	double tmp,df;
	double *res_ptr = res.dat;

	p2 = m2.dat;

	for (i2=0; i2< m2.ny; i2++)
	{
		p1 = m1.dat;
		for (i1=0; i1< m1.ny; i1++)
		{
			/* Inlining of dot */
			tmp =0;
			{
			  double *ap,*bp;
			  double *aend=&p1[m1.nx];
			  for (ap=&p1[0],bp=&p2[0]; ap<aend; ap++, bp++)
				  tmp+=*ap*(*bp);
			}
			(*res_ptr) = tmp;
			res_ptr++;
			p1+= m1.nx;
		}
		/* Advance to next column */
		p2 += m1.nx;
	}
}

void mat_el_exp(double *A, int nel)
{
	double *endp = &A[nel];
	double *p;

	for (p=A; p<endp; p++)
		(*p)=exp(*p);
}

double mat_el_exp_sum2(double *A, int nel)
{
	double *endp = &A[nel];
	double *p;
	double sum=0;

	for (p=A; p<endp; p++)
	{
		(*p)=exp(*p);
		sum += (*p);
	}

	return sum;
}


double mat_el_exp_sum(double *A, int nel)
{
	double *endp = &A[nel];
	double *p;
	double sum=0;
	double mx = -1e6;

	for (p=A; p<endp; p++)
	  if ((*p)>mx)
             mx=(*p);
	
	for (p=A; p<endp; p++)
	{
		(*p)=exp(*p-mx);
		sum += (*p);
	}
	return sum;
}

void mat_el_log(double *A, int nel)
{
	double ret=0;
	double *endp = &A[nel];
	double *p;

	for (p=A; p<endp; p++)
		(*p)=log(*p);
}


/* Add B into A */
void mat_add(double *A, double *B,int nel)
{
	double *endp = &A[nel];
	double *pa,*pb;

	for (pa=A,pb=B; pa<endp; pa++,pb++)
		(*pa)+=(*pb);
}

void mat_subtract(double *A, double *B,int nel)
{
	double *endp = &A[nel];
	double *pa,*pb;

	for (pa=A,pb=B; pa<endp; pa++,pb++)
		(*pa)-=(*pb);
}

void mat_el_prod(double *A, double *B,int nel)
{
	double *endp = &A[nel];
	double *pa,*pb;

	for (pa=A,pb=B; pa<endp; pa++,pb++)
		(*pa)*=(*pb);
}

void mat_el_mul_scalar(double *A, double s,int nel)
{
	double *endp = &A[nel];
	double *pa;

	for (pa=A; pa<endp; pa++)
		(*pa)*=s;
}

void mat_el_backdiv(double *A, double *B,int nel)
{
	double *endp = &A[nel];
	double *pa,*pb;

	for (pa=A,pb=B; pa<endp; pa++,pb++)
	{
		if ((*pb)==0 && (*pa)==0)
			(*pa)=1;
		else
			(*pa)=(*pb)/(*pa);
	}
}

void make_cond_dist(double *p,int n_dists,int var_size)

{
	double *prob_ptr,*dist_end;
	double dist_sum;
	int dist_i;

	prob_ptr = p;
	for (dist_i=0; dist_i<n_dists; dist_i++)
	{
		dist_end = prob_ptr+var_size;
		dist_sum=0;
		for (p=prob_ptr; p<dist_end; p++)
			dist_sum+= (*p);

		for (p=prob_ptr; p<dist_end; p++)
			(*p)/=dist_sum;
		prob_ptr+=var_size;
	}
}


double mat_max_abs(double *A, int nel)
{
    double curr;
    double max=-1e9;
	double *endp = &A[nel];
	double *p;

	for (p=A; p<endp; p++)
	{
		curr = (*p);
	    if (curr<0)
		    curr=-curr;
		if (curr>max)
    		max = curr;
		
    }  
	return max;
}


void write_mat(Matrix p)
{
	int ri,ci,ind=0;
	
    if (!bDbg)
        return;
	for (ri=0; ri<p.nx; ri++)
	{
		for (ci=0; ci<p.ny; ci++)
			printf("%g  ",GetMatEl(p,ri,ci));
		printf("\n");
	}
}

void write_vec(double *p,int n)
{
	int ri,ci,ind=0;
	
    if (!bDbg)
        return;
	for (ri=0; ri<n; ri++)
		printf("%g  ",p[ri]);
	printf("\n");
}

void dbgprintf(char *s)
{
    if (bDbg)
       printf(s);
}


/* This is actually minus the distance */
/* Generate the logemb model. a and b are vectors of length m1.ny and m2.ny respectively */

double get_maxel(Matrix phi, Matrix psi, double *a, double *b, int b_cond)

{
  int i1,i2,k,di;
  double *p1,*p2,*res_ptr;
  double tmp,df;
  double *ap,*bp;  
  double *aend;
  double logZ = 0;
  double max_el = -1e9;
  
  p2 = psi.dat;
  
  /* Reset everything */

  for (i2=0; i2< psi.ny; i2++)
  {
    p1 = phi.dat;
    for (i1=0; i1< phi.ny; i1++)
    {
      /* Inlining of l2_dist */
      tmp=0;
      {
	double d=0;
	aend=&p1[phi.nx];
	
	for (ap=&p1[0],bp=&p2[0]; ap<aend; ap++, bp++)	  {
	  d = *ap-*bp;
	  tmp+=(d*d);
	}
      }
      /* Calculate the factor in the exponent */
      tmp = -tmp+a[i1]+b[i2];
      if (tmp>max_el)
	max_el=tmp;
      logZ+= exp(tmp);
      p1+= phi.nx;
    }
    p2+= phi.nx;
  }
  logZ = log(logZ);
  /*  return logZ; */
  return max_el;
}

/* This is actually minus the distance */
/* Generate the logemb model. a and b are vectors of length m1.ny and m2.ny respectively */

double get_worst_max(Matrix phi, Matrix psi, double *a, double *b, int b_cond)
{
  double max_el;
  double phi_max,psi_max,a_max,b_max;
  
  phi_max =  mat_max_abs(phi.dat, phi.nall);
  psi_max =  mat_max_abs(psi.dat, psi.nall);

  /* Calculate upper bound on distance between phi and psi */
  if (phi_max>psi_max)
    max_el = 4*phi_max*phi_max*phi.nx;
  else
    max_el = 4*psi_max*psi_max*phi.nx;

  a_max = mat_max_abs(a, phi.ny);
  b_max = mat_max_abs(b, psi.ny);  
  if (a_max>b_max)
    max_el+=a_max;
  else
    max_el+=b_max;

  return max_el;
}

/* This is actually minus the distance */
/* Generate the logemb model. a and b are vectors of length m1.ny and m2.ny respectively */

void mat_logemb_model(Matrix phi, Matrix psi, double *a, double *b, Matrix phi_expect, Matrix psi_expect,  double *px, double *py
		      ,double *p0_i, double *p0_j, double *p0_val, int n_p0, double *loglik, int b_cond)

{
  int i1,i2,k,di;
  double *p1,*p2,*res_ptr;
  double *p_phiexp,*p_psiexp;
  double tmp,df;
  double Z=0;
  int D = phi_expect.nx;
  double *ap,*bp,*p_phie,*p_psie;  
  int NX = phi.ny;
  int NY = psi.ny;
  int p0ind = 0;
  double *aend;
  double logZ;
  
  (*loglik) = 0;
  p2 = psi.dat;
  p_phiexp = phi_expect.dat;
  p_psiexp = psi_expect.dat;
  
  /* Reset everything */

  /* Initialize max element only if needed */
  if (max_element==0)
    /*    max_element = get_worst_max(phi, psi, a, b, b_cond);    */
    max_element = get_maxel(phi, psi, a, b, b_cond);
 
  if (!b_cond)
  {
    memset(phi_expect.dat,0,phi_expect.nall*sizeof(double));
    memset(psi_expect.dat,0,psi_expect.nall*sizeof(double));
  }
  memset(px,0,NX*sizeof(double));
  memset(py,0,NY*sizeof(double));

  
  for (i2=0; i2< psi.ny; i2++)    {
    p1 = phi.dat;
    p_psiexp = psi_expect.dat;
    
    for (i1=0; i1< phi.ny; i1++)     {
	double d=0;
        
      /* Inlining of l2_dist */
      tmp=0;
      {
	aend=&p1[phi.nx];
	
	for (ap=&p1[0],bp=&p2[0]; ap<aend; ap++, bp++)	  {
	  d = *ap-*bp;
	  tmp+=(d*d);
	}
      }
      tmp = -tmp+a[i1]+b[i2]-max_element;

      if ((i1==p0_i[p0ind]-1) && (i2==p0_j[p0ind]-1))
      {
	(*loglik)+= p0_val[p0ind]*tmp;
	p0ind++;
      }
      
      tmp = exp(tmp);
      
      if (!b_cond)
      {
	/* This can be speeded up by avoiding called to GetMatInd */
	aend=&p1[D];
	for (di=0,ap=&p1[0],bp=&p2[0],p_phie=&p_phiexp[0],p_psie=&p_psiexp[0]; di<D; di++,ap++,bp++,p_phie++,p_psie++)
        {
	  /*
	    phi_expect.dat[GetMatInd(phi_expect,di,i2)]+=tmp*(*ap);
	    psi_expect.dat[GetMatInd(psi_expect,i1,di)]+=tmp*(*bp);
	  */
	  /*	
	    phi_expect.dat[GetMatInd(phi_expect,di,i2)]+=tmp*phi.dat[GetMatInd(phi,di,i1)];
	    psi_expect.dat[GetMatInd(psi_expect,i1,di)]+=tmp*psi.dat[GetMatInd(psi,di,i2)];
	  */
	  (*p_phie)+=tmp*(*ap);
	  (*p_psie)+=tmp*(*bp);
	}
      }
      px[i1]+=tmp;
      py[i2]+=tmp;

      Z+= tmp;

      p1+= phi.nx;
      p_psiexp += D;      
    }
    /* Advance to next column */
    p2 += phi.nx;
    p_phiexp += D;          
  }

  /* Verify normalization to one, since we divided by the true logZ */
  if (!b_cond)
  {
     vec_mul_scalar(px,NX,1/Z);
     vec_mul_scalar(py,NY,1/Z);
     mat_el_mul_scalar(phi_expect.dat, 1/Z, phi_expect.nall);
     mat_el_mul_scalar(psi_expect.dat, 1/Z, psi_expect.nall);
     (*loglik)-= log(Z);
  }
}


void AllocMat(int nx, int ny, Matrix *m)
{
  m->nx = nx;
  m->ny = ny;
  m->nall = nx*ny;
  m->dat = malloc(sizeof(double)*nx*ny);
  if (m->dat==NULL)
    printf("Error in allocation\n");

}

void mat_sum_rows(Matrix m, double *sum)
{
	int xi,yi;
	for (xi=0; xi<m.nx; xi++)
	{
		sum[xi]=0;
		for (yi=0; yi<m.ny; yi++)
			sum[xi]+=GetMatEl(m,xi,yi);

	}
}

#if 0 
void mat_sum_cols(Matrix m, double *sum)
{
	int xi,yi;

	for (yi=0; yi<m.ny; yi++)
	{
		sum[yi]=0;
		for (xi=0; xi<m.nx; xi++)
			sum[yi]+=GetMatEl(m,xi,yi);
	}
}
#endif


void mat_sum_cols(Matrix m, double *sum)
{
	int xi,yi;
	double *p=m.dat;

	for (yi=0; yi<m.ny; yi++,p+=m.nx)
	{
		sum[yi]=vec_sum(p,m.nx);
	}
}

/* Multiply columns of m by v */
void mat_mul_cols(Matrix m, double *v)
{
	int xi,yi;
	double *p=m.dat;

	for (yi=0; yi<m.ny; yi++,p+=m.nx)
		vec_mul_scalar(p,m.nx,v[yi]);
}

void mat_mul_rows(Matrix m, double *v)
{
	int xi,yi;
	for (xi=0; xi<m.nx; xi++)
	{
		for (yi=0; yi<m.ny; yi++)
			m.dat[GetMatInd(m,xi,yi)]*=v[xi];

	}
}

void MatFromMLab(const mxArray *plhs, Matrix *mat)
{
	mat->dat = mxGetPr(plhs);
	mat->nx = mxGetM(plhs);
	mat->ny = mxGetN(plhs);
	mat->nall = mat->nx*mat->ny;
}

void mat_set(Matrix A, Matrix B)
{
  memcpy(A.dat,B.dat,A.nall*sizeof(double));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double M;
  double log_thr;
  int Xs,D,Ys,i,j,ep=0;
  int maxentep;
  int pxy_size;
  double mx;
  Matrix phi,psi,phi_tr,psi_tr,phi_exp_d,psi_exp_d,psi_exp_m,phi_exp_m,grad_phi,grad_psi;
  double *px0,*py0,*px,*py,*pa,*pb,*grad_a,*grad_b;
  double Z,loglik;
  int xi,yi,di;
  int b_cond; /* Says if the model is conditioned on x, or joint */
  double *p0_i,*p0_j,*p0_val;
  int Np0;
  
  MatFromMLab(prhs[0],&phi);   /* Size NX,D */ 
  MatFromMLab(prhs[1],&psi);   /* Size NX,D */	
  MatFromMLab(prhs[2],&phi_tr);   /* Size D,NX */ 
  MatFromMLab(prhs[3],&psi_tr);   /* Size D,NY */	
  px0		  = mxGetPr(prhs[4]);
  py0		  = mxGetPr(prhs[5]);
  MatFromMLab(prhs[6],&phi_exp_d);
  MatFromMLab(prhs[7],&psi_exp_d);
  
  pa = mxGetPr(prhs[8]);
  pb = mxGetPr(prhs[9]);
  b_cond = *mxGetPr(prhs[10]);

  p0_i = mxGetPr(prhs[11]);
  p0_j = mxGetPr(prhs[12]);
  p0_val = mxGetPr(prhs[13]);
  
  
  Xs = mxGetM(prhs[0]);
  D = mxGetN(prhs[0]);
  Ys = mxGetM(prhs[1]);
  Np0 = mxGetNumberOfElements(prhs[13]);
    
  pxy_size = Xs*Ys;
  
  plhs[0] = mxCreateDoubleMatrix(D,Xs,mxREAL); 
  plhs[1] = mxCreateDoubleMatrix(D,Ys,mxREAL);
  plhs[2] =  mxCreateDoubleMatrix(1,Xs,mxREAL);
  plhs[3] =  mxCreateDoubleMatrix(1,Ys,mxREAL);
  plhs[4] =  mxCreateDoubleMatrix(1,1,mxREAL);
  
  grad_phi.dat = mxGetPr(plhs[0]);
  grad_psi.dat = mxGetPr(plhs[1]);
  grad_phi.nx = D;
  grad_phi.ny = Xs;
  grad_phi.nall = Xs*D;
  
  grad_psi.nx = D;
  grad_psi.ny = Ys;
  grad_psi.nall = Ys*D;
  
 
  grad_a = mxGetPr(plhs[2]);
  grad_b = mxGetPr(plhs[3]);
  
  
  AllocMat(D,Ys,&phi_exp_m);
  AllocMat(D,Xs,&psi_exp_m);
  
  px  = malloc(Xs*sizeof(double)); 
  py  = malloc(Ys*sizeof(double)); 

  if (px==NULL || py==NULL)
    printf("Error in allocation\n");

  max_element = 0;
  
  if (!b_cond)
    mat_logemb_model(phi_tr, psi_tr, pa, pb, phi_exp_m, psi_exp_m, px, py, p0_i, p0_j, p0_val, Np0, mxGetPr(plhs[4]),b_cond) ;
  else
  {
    /* Run logemb model one time to get the x partition function in px */
    mat_logemb_model(phi_tr, psi_tr, pa, pb, phi_exp_m, psi_exp_m, px, py, p0_i, p0_j, p0_val, Np0, mxGetPr(plhs[4]),b_cond) ;
    /* Create new x factor */
    for (xi=0; xi<Xs; xi++)
        grad_a[xi] =   pa[xi]-log(px[xi])+log(px0[xi]);
    mat_logemb_model(phi_tr, psi_tr, grad_a, pb, phi_exp_m, psi_exp_m, px, py, p0_i, p0_j, p0_val, Np0, mxGetPr(plhs[4]),0) ;
  }


  /* Matlab line: diag(px)*phi-pxy*psi - (diag(px0)*phi-pxy0*psi);
     /* equivalent to  diag(px-px0)*phi+psi_expect_d - psi_expect_m */
  vec_subtract(px,px0,Xs);
  mat_set(grad_phi,phi_tr);
  /*  printf("\n Phi1 Max=%g MaxPx=%g\n",mat_max_abs(grad_phi.dat,grad_phi.nall),mat_max_abs(px,Xs));    */
  mat_mul_cols(grad_phi,px);
  /*  printf("\n Phi2 Max=%g MaxPx=%g\n",mat_max_abs(grad_phi.dat,grad_phi.nall),mat_max_abs(px,Xs)); */
  mat_subtract(psi_exp_d.dat,psi_exp_m.dat,psi_exp_d.nall);
  /* printf("\n Phi3 Max=%g\n",mat_max_abs(grad_phi.dat,grad_phi.nall));      */
  mat_add(grad_phi.dat,psi_exp_d.dat,psi_exp_d.nall);
  /* printf("\n Phi4 Max=%g\n",mat_max_abs(grad_phi.dat,grad_phi.nall)); */
  
  vec_subtract(py,py0,Ys);
  mat_set(grad_psi,psi_tr);
  mat_mul_cols(grad_psi,py);
  mat_subtract(phi_exp_d.dat,phi_exp_m.dat,phi_exp_d.nall);
  mat_add(grad_psi.dat,phi_exp_d.dat,phi_exp_d.nall);

  /* Generate gradient for a and b. Remember px already has px  subtracted. Just need to flip sign */
  memcpy(grad_a,px,Xs*sizeof(double));
  vec_mul_scalar(grad_a,-1,Xs);
  
  memcpy(grad_b,py,Ys*sizeof(double));
  vec_mul_scalar(grad_b,-1,Ys);
  
  free(psi_exp_m.dat);
  free(phi_exp_m.dat);
  free(px);
  free(py);
  
}


