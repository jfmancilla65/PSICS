# psics python package
# file: call_glasso.R
# R function to call glasso algorithm

# R glasso library 
library(glasso)

# 1) read parameters from command line

args = commandArgs(TRUE)

# 2) glasso execute parameters read from command line
sigma_file = args[1]
rho        = args[2]
nobs       = args[3] 
zero       = args[4] 
thr        = as.double(args[5])
maxit      = as.numeric(args[6])
approx     = args[7]
penalize.diagonal = args[8]
start      = args[9]
w.init     = args[10]
wi.init    = args[11]
trace      = args[12]


# 3) read covariance matrix data and rho value(s)
sigma = read.csv(sigma_file, header = FALSE)
sigma = as.matrix(sigma)

if (rho == "_glassorho.dat") {
 rho = read.csv("_glassorho.dat", header = FALSE)
 m = nrow(rho)
 n = ncol(rho)
 if (n == 1) {
  rho = as.vector(unlist(rho))
 } else {
  rho = as.matrix(rho)
 }
} else {
   rho = as.numeric(rho)
}

# 4) input arguments processing
if (trace == "TRUE"){
 trace = TRUE
} else {
   trace = FALSE
}

if (nobs == "NULL") {
 nobs = NULL
} else { 
   nobs = as.numeric(nobs)
}

if (zero == "NULL") { 
 zero = NULL
} else {
	 zero = read.csv(zero, header = FALSE)
	 zero = as.matrix(zero)
}

if (approx == "FALSE") {
 approx = FALSE
} else {
   approx = TRUE   
}

if (penalize.diagonal == "FALSE") {
 penalize.diagonal = FALSE
} else {
   penalize.diagonal = TRUE   
}

if (start == "warm") {
 # load initial matrices files
 w.init  = read.csv(w.init,  header = FALSE)
 w.init  = as.matrix(w.init)
 wi.init = read.csv(wi.init, header = FALSE)
 wi.init = as.matrix(wi.init) 
}

# 5) call glasso function
gs = glasso(sigma,rho,nobs=nobs,zero=zero,thr=thr,maxit=maxit,approx=approx,
            penalize.diagonal=penalize.diagonal,start=start,w.init=w.init,wi.init=wi.init,trace=trace)

# 6) save results in temp files 
XSol    = gs$wi
USol    = gs$w
loglik  = gs$loglik
errflag = gs$errflag
approx  = gs$approx
del     = gs$del
niter   = gs$niter
result = list(c('loglik','errflag','approx','del','niter'),
              c(loglik,errflag,approx,del,niter))
write.table(XSol, file ="_glassoXsol.dat", sep = ",", col.names = FALSE,row.names=FALSE)
write.table(USol, file ="_glassoUsol.dat", sep = ",", col.names = FALSE,row.names=FALSE)
write.table(result, file ="_glassores.dat",  sep = ",", col.names = FALSE,row.names=FALSE,quote=FALSE)


