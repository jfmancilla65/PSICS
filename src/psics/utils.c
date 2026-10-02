#include "BoxQP.h"

// BLAS wrapper for win32 version
//#ifdef win32
/*
void cblas_dscal( int N, double alpha, double *X, int incX)
{
	dscal(&N,&alpha,X,&incX);
}

void cblas_dcopy(int N,double *X,int incX,double *Y,int incY)
{
	dcopy(&N,X,&incX,Y,&incY);
}

void cblas_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB,
                 int M, int N, int K, double alpha, double *A, int lda,
                 double *B, int ldb, double beta, double *C, int ldc)
{
	char ta[1],tb[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	if (transB==111)
	{
		*tb='N';
	}
	else
	{
		*tb='T';
	};
	dgemm(ta,tb,&M,&N,&K,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

void cblas_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA,
                 int M, int N, double alpha, double *A, int lda,
                 double *B, int incB, double beta, double *C, int incC)
{
	char ta[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	dgemv(ta,&M,&N,&alpha,A,&lda,B,&incB,&beta,C,&incC);
}

void cblas_daxpy(int N,double alpha,double *X,int incX,double *Y,int incY)
{
	daxpy(&N,&alpha,X,&incX,Y,&incY);
}

//#endif win32

*/

// Some useful functions ...
double doubsum(double *xmat, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xmat[i];};
	return res;
}


double doubdot(double *xvec, double *yvec, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xvec[i]*yvec[i];};
	return res;
}

int idxmax(double *xmat, int n)
{
	int i;
	int res=0;
	for (i=0;i<n;i++)
	{
		if (xmat[i]>xmat[res]) {res=i;}
	}
	return res;
}


double doubasum(double *xmat, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=dabsf(xmat[i]);};
	return res;
}

double doubnorm2(double *xmat, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xmat[i]*xmat[i];};
	return sqrt(res);
}

double dsignf(double x)
{
	if (x>=0)
		return 1.0;
	else
		return -1.0;
}

double dminif(double x, double y)
{
	if (x>=y)
		return y;
	else
		return x;
}

int imaxf(int x, int y)
{
	if (x>=y)
		return x;
	else
		return y;
}

double dabsf(double x)
{
	if (x>=0)
		return x;
	else
		return -x;
}

void dispmat(double *xmat, int n, int m)
{
	int i,j;

	for (i=0; i<n; i++)
	{
		for (j=0;j<m;j++)
		{
			printf("%+.4f ",xmat[j*n+i]);
		}
		printf("\n");
	}
	printf("\n");
}

double maxeig(double *xmat, int n) 
{
	// xmat is symmetric n x n matrix
	int incx=1,indmax,maxloop=10000,k=0;
	double alpha, beta, dmax, dmax_temp,dmax_tol;
	double *bufveca=(double *) calloc(n,sizeof(double));
	double *bufvecb=(double *) calloc(n,sizeof(double));

	dmax_tol=.001;
	// do power series to get approximation to max eigenvalue of A+X
	alpha=0.0;cblas_dscal(n,alpha,bufveca,incx);bufveca[0]=1.0; // x_0 = [1,0,0,...] do something better later
	beta=0.0;dmax=1.0;dmax_temp=0.0;
	while ((dabsf(dmax-dmax_temp)>dmax_tol)&&(k<=maxloop)){
		dmax_temp=dmax;
		alpha=1.0;
		cblas_dgemv(CblasColMajor,CblasNoTrans,n,n,alpha,xmat,n,bufveca,incx,beta,bufvecb,incx);
		indmax=idxmax(bufvecb,n);dmax=bufvecb[indmax];
		alpha=1.0/dmax;cblas_dscal(n,alpha,bufvecb,incx);
		cblas_dcopy(n,bufvecb,incx,bufveca,incx);
		k++;
	}
	alpha=1.0;	
	// compute Rayleigh Quotient to approximate max eigenvalue of A+X
	cblas_dgemv(CblasColMajor,CblasNoTrans,n,n,alpha,xmat,n,bufvecb,incx,beta,bufveca,incx);
	dmax=doubdot(bufvecb,bufveca,n);alpha=doubnorm2(bufvecb,n);dmax=dmax/alpha/alpha;
	free(bufveca);free(bufvecb);
	return dmax;
}
