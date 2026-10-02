# PSICS

## A Python toolbox for numerically solving the sparse inverse covariance selection. 

## Introduction

This Python module includes integrated access to native code for **COVSEL**, **COVPATH**, **GLASSO**, and **ALM** algorithms that were previously published. 

## PSICS description

PSICS is built using the Python programming language.  Consists of several functions that integrate the call to each of the sparse inverse covariance selection (SICS) algorithms, using the standard library **subprocess**.

This approach makes it possible to call different software scripts or programs written in **Matlab** and **R**, the two main platforms used to build original algorithms.

PSCIS installs the native code algorithms **COVSEL**, **COVPATH**, **GLASSO** and **ALM**. Originally, **COVSEL** and ALM runned in a **Matlab** environment. In this package, **Matlab** was replaced with the Octave platform.  

Detailed directions to install the PSICS package can be found in the README.md file. 

## PSICS functions

The following sections describe the python functions that call the SICS algorithms, including input and output parameters.

### COVSEL

``` 
# python function call
from psics import covsel
outcovsel = covsel.call_covsel(sigma,rho=0.01,maxiter=100,prec=1.e1,maxnest=100,algo='BoxQP')
```

**Input parameters**
|Parameter     | Description| Type      | Required | Default |
|--------------|-----------|:----------:|:--------:|:-------:|
|sigma|Covariance matrix (symmetric)| array(n,n) float|y||
|rho|Sparsity regularization factor| float|y|0.5|
|maxiter|Number of iterations in main loop|int|y|100|
|prec|Target precision|float|y|1e-1|
|maxnest|Maximum number of iterations in quadratic program solver| int|y|100|
|algo|Quadratic program solver to use|string, ("BoxQP","Sedumi")|y|"BoxQP"|

**Output dictionary parameters**

|Parameter     | Description| Type      |
|--------------|-----------|:----------:| 
|outcovsel.X|Estimated inverse covariance matrix  |array(n,n), float| 
|outcovsel.U|Optimal positive-definite matrix| array(n,n), float|
|outcovsel.gvals|Gap values in dual per iteration| array(m), float|
|outcovsel.cputimes|CPU time per iteration in seconds| array(m), float|

### GLASSO


The code for this algorithm is an R package available in the CRAN repository to download or install directly in the user R environment. It requires a FORTRAN compiler in order to link code to the package. All processes to set up the R GLASSO package are done automatically by the PSICS installer. The main R function that implements this algorithm is glasso(). The Python function shown next is a wrapper of the R function call.



```
# python function call
from psics import glasso

outglasso = glasso.call_glasso(sigma, rho=0.001, nobs=NULL, zero=NULL, thr=1.0e-4, maxit=1e4,penalizediagonal=TRUE, start="cold", winit=NULL,wiinit=NULL, trace=FALSE)

```

**Input parameters**

|Parameter     | Description| Type      | Required | Default |
|--------------|-----------|:----------:|:--------:|:-------:|
|sigma | Covariance  matrix (symmetric)| array(n,n), float| y | |
| rho | Non-negative regularization parameter for lasso. rho=0 means no regularization. Can be a scalar (usual)  or a vector of length n | float or array(n)|y |0.5| 
| thr          |Threshold for convergence. Iterations stop when average absolute parameter change is less than thr*ave(abs(offdiag(s))) | int      | y          | 1e-4         |
|start          |Using start="warm"  can provide starting values for w and wi| string, ("cold","warm")|y| "cold"| 
|rho | Non-negative regularization parameter for lasso. rho=0 means no regularization. Can be a scalar (usual)  or a vector of length|float or array(n)|y|0.01|
|nobs|Number of observations used in computation of the covariance matrix s. This quantity is need to compute the value of log-likelihood. If not specified, loglik will be returned as Null.|int|n|NULL|
|zero|Indices of entries of inverse covariance to be constrained to be zero. The input should be a matrix with two columns, each row indicating the indices of elements to be constrained to be zero. The solution must be symmetric, so you need only specify one of (j,k) and (k,j). An entry in the zero matrix overrides any entry in the rho matrix for a given element.|array(j,k) float|n|NULL|
|thr|Threshold for convergence. Iterations stop when average absolute parameter change is less than thr*ave(abs(offdiag(s)))| float|y|1e-4|
|maxit|Maximum number of iterations of outer loop. | int|y| 10,000| 
|approx|  Approximation flag: if TRUE, computes Meinhausen-Buhlmann(2006) approximation |Boolean|y|FALSE|
|penalizediagonal|Should diagonal of inverse covariance be penalized?|Boolean|y|TRUE|
|start| Can provide starting values for w and wi |String ("cold","warm")|y|"cold"|
|winit|Optional starting values for estimated covariance matrix. Only needed when start="warm" is specified|array(n,n)|n|NULL| 
|wiinit|Optional starting values for estimated inverse covariance. Only needed when start="warm" is specified| array(n,n)|n|NULL| 
|trace|Flag for printing out information as iterations proceed|boolean|y|FALSE


**Output dictionary parameters**

|Parameter     | Description| Type      |
|--------------|-----------|:----------:| 
|outglasso.wi|Estimated  inverse covariance matrix|array(n,n) float|   
|outglasso.w| array(n,n) float|Estimated  covariance matrix|
|outglasso.loglik|Value of maximized log-likelihood+penalty|float|
|outglasso.errflag|Error flag. Zero value if no memory allocation errors ocurred|int|
|outglasso.approx|Value of input argument approx| boolean|
|outglasso.del|Change in parameter value at convergence|float|
|outglasso.niter|Number of iterations of outer loop used by algorithm|int|



### GLASSOPATH

This is a variant of the GLASSO algorithm and implements a pathwise approach to solve the SICS problem. Its executable code is included in the same R package, and the main function that executes the algorithm is glassopath(). The wrapper Python function that invokes the R code is shown next.

```
# python function call
from psics import glasso

outglassopath =  glasso.call_glassopath(sigma,rholist=rholist,trace=0,thr=1e-, maxit=1e4, approx=FALSE, penalizediagonal=TRUE, winit=NULL,wiinit=NULL, trace=1)
```

**Input parameters**

|Parameter     | Description| Type      | Required | Default |
|--------------|-----------|:----------:|:--------:|:-------:|
|sigma|Covariance positive-definite matrix (symmetric)|array(n,n), float|y|| 
|rholist|Regularization parameters for the lasso. Should be increasing from smallest to largest; actual path is computed from largest to smallest value of rho. If NULL, 10 values in a (hopefully reasonable) range are used. Note that the same parameter rholist[j] is used for all entries of the inverse covariance matrix; different penalties for different entries are not allowed.|array(m), float|y|NULL|
|thr| Threshold for convergence. Iterations stop when average absolute parameter change is less than thr*ave(abs(offdiag(s)))| float |y|1e-4|
|maxit|Maximum number of iterations of outer loop| int|y|10,000|
|approx|Approximation flag: if true, computes Meinhausen-Buhlmann(2006) approximation|boolean|y|FALSE|  
|penalizediagonal|Should diagonal of inverse covariance be penalized?|boolean|y|TRUE|
|winit|Starting values for estimated covariance matrix. Only needed when start="warm" is specified| array(n,n)|n|NULL|
|wiinit|Optional starting values for estimated inverse covariance matrix. Only needed when start="warm" is specified| array(n,n)|n|NULL|
|trace| Flag for printing out information as iterations proceed.  If trace=0 means no printing, trace=1 means outer level printing, trace=2 means full printing| int, (0,1,2)|y|0|

**Output dictionary parameters**

|Parameter     | Description| Type      |
|--------------|-----------|:----------:| 
|outglassopath.wi|Estimated inverse covariance matrices|array(n,n,size(rholist)), float.|
|outglassopath.w|Estimated  covariance matrices|array(n,nsize(rholist)), float|
|outglassopath.errflag|Error flag, zero value if no errors ocurred|int|
|outglassopath.approx|Value of input argument approx|float| 
|outglassopath.rholist|Values of regularization parameter used| array(p)|
|outglassopath.errflag|Values of error flag, 0 means no memory allocation error|int| 



### COVPATH

COVPATH is also coded in the R platform and is necessary to link some C functions. The main R function to execute this algorithm is covpath() and the Python routine that invokes this function is described here.

```
# python function call
from psics import covpath
outcovpath  = covpath.call_covpath(sigma,rholist=rholist,t=1e-2,colFraction=1)
```


**Input parameters**

|Parameter     | Description| Type      | Required | Default |
|--------------|-----------|:----------:|:--------:|:-------:|
|sigmaSample covariance positive-definite matrix|, array(n,n), float. |y||
|rholist|Sparsity regularization factor, ordered from largest to smallest value of rho. If NULL, 10 values in a (hopefully reasonable) range are used.|array(m), float|y|NULL| 
|t|Factor that controls smoothness| float|y|1e-7|
|colFraction|Blocksize used in coordinate descent algorithm|int|y|1| 



**Output dictionary parameters**

|Parameter     | Description| Type      |
|--------------|-----------|:----------:| 
|outcovpath.invCovmat|Estimate of inverse covariance matrix|array(n,n), float|  
|outcovpath.Covmat|Estimate of covariance matrix|array(n,n), float| 
|outcovpath.relativeGaps|Relative gaps|array(n), float| 



**ALM**

The ALM algorithm code is written in Matlab and the main function is found in file SICS_ALM.m. This code is called using the following toolbox Python function.

```
from psics import alm
outalm = alm.call_alm(sigma,rho=0.5,mxitr=100,mu0=2.0,muf=0.0001,rmu=0.25,tol_gap=0.1,record=1,fsigma=1e2)
```

**Input parameters**

|Parameter     | Description| Type      | Required | Default |
|--------------|-----------|:----------:|:--------:|:-------:|
|sigma|Covariance matrix|array(n,n), float|y| 
|rho|Sparsity regularization factor|float|y|0.5|
|mxitr|Maximum number of iterations|int|y|500| 
|mu0|Mu initial value|float|y|1e-1| 
|muf|Mu final value|float|y|1e-3|
|rmu|Mu ratio of decreasing mu|float|y|0.25| 
|tol_gap|Tolerance for duality gap|float|y|1e-1| 
|tol_frel|Tolerance for relative change of obj value|float|y|1e-7|
|tol_Xrel|Tolerance for relative change of X|float|y|1e-7| 
|tol_Yrel|Tolerance for relative change of Y|float|y|1e-7| 
|tol_pinf|Tolerance for infeasibility|float|y|1e-3| 
|numDG|Every numDG iterations, we compute duality gap since it's expensive|int|y|10|
|record|Print out information if record > 0|int|y|1| 
|fsigma|Smoothness parameter| float|y|1e-10|


**Output dictionary parameters**

|Parameter     | Description| Type      |
|--------------|-----------|:----------:| 
|alm.X|Estimated inverse covariance matrix|array(n,n), float| 
|alm.Y|Y matrix|array(n,n), float| 
|alm.itr|Number of iterations executed|int| 
|alm.obj|Objective function value|float|
|alm.gap|Dual gap|float| 
|alm.gapX|Dual gap X|float|
|alm.gapY|Dual gap Y|float|
|alm.pinf|P inf value|float|







