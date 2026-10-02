# psics python package
# file: call_covpath.R
# calls covpath algorithm 

library(Covpath)

args = commandArgs(TRUE)

# read data from command line 
sigma_file   = args[1]
rholist_file = args[2]
rhomax_s     = args[3]
t            = as.double(args[4])
colFraction  = as.numeric(args[5])


# covariance matrix
Sigma = read.csv(sigma_file, header = FALSE)
Sigma = as.matrix(Sigma)

# rhomax value
if(rhomax_s == "NULL") {
 rhomax = max(diag(Sigma))
} else {
 rhomax = as.numeric(rhomax_s)
}

# rholist default vector
if (rholist_file == "NULL") {
 nrhos = 10
 logrhomax = log10(rhomax)
 logrhomin = logrhomax - 1.5
 rholist = 10^(seq(from = logrhomax, to = logrhomin, length.out= (nrhos+1)))
} else {
 rholist = read.csv(rholist_file, header = FALSE)
 rholist = as.vector(unlist(rholist))
}

# sort rholist 
rholist = sort(rholist,decreasing = TRUE)
n  = nrow(Sigma)



# Covpath Parameters

covpathsol = covpath(Sigma,rholist,t,colFraction,rhomax)
XSol   = covpathsol$invCovmat
USol   = covpathsol$Covmat
RelGap = covpathsol$relativeGaps

# test, save all original matrix results
write.table(XSol,    file ="_covpathXsol.dat", sep = ",", col.names = FALSE,row.names=FALSE)
write.table(USol,    file ="_covpathUsol.dat", sep = ",", col.names = FALSE,row.names=FALSE)
write.table(RelGap,  file ="_covpathRelGap.dat", sep = ",", col.names = FALSE,row.names=FALSE)
write.table(rholist, file ="_covpathRholist.dat", sep = ",", col.names = FALSE,row.names=FALSE)
