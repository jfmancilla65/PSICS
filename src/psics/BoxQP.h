//#define win32
#ifdef MACVERSION 
#define mac 
#endif

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mex.h>
#include <cblas.h>


#ifdef win32
#include <malloc.h>
#endif

#ifdef mac
#include "/System/Library/Frameworks/vecLib.framework/Headers/cblas.h"
#include "/System/Library/Frameworks/vecLib.framework/Headers/clapack.h"
#endif

// Calling LAPACK on the mac
#ifdef mac
#define dsyev dsyev_ 
#define dgesvd dgesvd_
#endif

#ifdef win32
enum CBLAS_ORDER 	{CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE 	{CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO		{CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG		{CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE		{CblasLeft=141, CblasRight=142};
void cblas_dscal(int N,double alpha, double *X,int incX);
void cblas_dcopy(int N,double *X,int incX,double *Y,int incY);
void cblas_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);
void cblas_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, int M, int N, double alpha, double *A, int lda, double *B, int incB, double beta, double *C, int incC);
void cblas_daxpy(int N,double alpha,double *X,int incX,double *Y,int incY);
double cblas_ddot(int n, double *x, int incx, double *y, int incy); 
#endif


// Local functions
double doubsum(double *xmat, int n);

double doubdot(double *xvec, double *yvec, int n);

double doubasum(double *xmat, int n);

double doubnorm2(double *xmat, int n);

int idxmax(double *xmat, int n);

void BoxQP(double *Amat, int n, double *lvec, double *uvec, int MaxIter, int info, int verbose, double *xvecout);

double dsignf(double x);

double dminif(double x, double y);

int imaxf(int x, int y);

double dabsf(double x);

void dispmat(double *xmat, int n, int m);

double maxeig(double *xmat, int n); 
