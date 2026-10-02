# psics python package
# file:call_glassopath.R
# R function to call glassopath algorithm

# R glasso library 
library(glasso)

args = commandArgs(TRUE)

# 1) read data from command line
sigma_file = args[1]
sigma = read.csv(sigma_file, header = FALSE)
sigma = as.matrix(sigma)
n  = nrow(sigma)
rholist = args[2]
thr     = as.double(args[3])
maxit   = as.numeric(args[4])
approx  = args[5]
penalize.diagonal = args[6]
w.init  = args[7]
wi.init = args[8]
trace   = as.numeric(args[9])


# 2) input arguments processing

# rholist 
if (rholist != "NULL") {
 rholist = read.csv(rholist, header = FALSE)
 rholist = as.vector(unlist(rholist))
 rholist = sort(rholist,decreasing = FALSE)
} else { rholist = NULL }

# approx
if (approx == "TRUE")  { 
 approx = TRUE 
} else { approx = FALSE }

#penalized diagonal
if (penalize.diagonal == "FALSE") {
 penalize.diagonal = FALSE
} else { penalize.diagonal = TRUE }

# w.init 
if (w.init != "NULL") {
 w.init = read.csv(w.init, header = FALSE)
 w.init  = as.matrix(w.init)
} else { w.init = NULL }

# wi.init
if (wi.init != "NULL") {
 wi.init = read.csv(wi.init, header = FALSE)
	wi.init = as.matrix(wi.init)  
} else { wi.init = NULL}

# 3) call glasso function

gs = glassopath(sigma,rholist=rholist,thr=thr,maxit=maxit,approx=approx,
                penalize.diagonal=penalize.diagonal,w.init=w.init,wi.init=wi.init,trace=trace)

# 4) save results in temp files 
wi      = gs$wi
w       = gs$w
errflag = gs$errflag
approx  = gs$approx
rholist = gs$rholist


write.table(rholist, file ="_glassopathrholist.dat", sep = ",", col.names = FALSE,row.names=FALSE)
write.table(errflag, file ="_glassopatherrflag.dat", sep = ",", col.names = FALSE,row.names=FALSE)
results = list(c('approx'),c(approx))
write.table(results, file ="_glassopathres.dat",  sep = ",", col.names = FALSE,row.names=FALSE,quote=FALSE)


# resulting matrices are stored in a row for each value of rholist 
m = length(rholist)
X = matrix(0,nrow=m,ncol=n*n)
S = matrix(0,nrow=m,ncol=n*n)
 
for (ind in(m:1))
{
	k = 1
  for (i in 1:n) {
   for (j in 1:n) {
      X[ind,k] = wi[i,j,ind]  # inverse of covariance matriz for rho[ind]
      S[ind,k] = w[i,j,ind]   # covariance matrix for rho[ind]
      k = k +1
   }
  }
}

write.table(X, file ="_glassopathXsol.dat", sep = ",", col.names = FALSE,row.names=FALSE)
write.table(S, file ="_glassopathUsol.dat", sep = ",", col.names = FALSE,row.names=FALSE)

