# psics python package
# file: utils.py
# support functions
#
import numpy as np
import math
import numpy.linalg as LA

#################################################################################### 
def d1(A,B):
# simple distance between matrices A, B 
####################################################################################
 ma, na = A.shape
 mb, nb = B.shape
 dv = 0
 if ma == mb and na == nb:
  for i in range(ma):
   for j in range(i):
    dv += np.abs(A[i,j]-B[i,j])
 else:
  dv = -1
 return dv
 
#################################################################################### 
def d2(A,B):
# Euclidian distance between matrices A, B
####################################################################################
 ma, na = A.shape
 mb, nb = B.shape
 dv = 0
 if ma == mb and na == nb:
  for i in range(ma):
   for j in range(i):
    dv += (A[i,j] - B[i,j]) * (A[i,j] - B[i,j])
  dv = math.sqrt(dv) 
 else:
  dv = -1
 return dv
 
#################################################################################### 
def d3(A,B):
####################################################################################
 ma, na = A.shape
 mb, nb = B.shape
 dv = 0
 if ma == mb and na == nb:
  logA = np.log(np.sum(A))
  logB = np.log(np.sum(B))
  for i in range(ma):
   for j in range(i):
    dv += A[i,j] * ( np.log(A[i,j]/B[i,j]) + logB - logA)
 else:
  dv = -1
 return dv
 
####################################################################################   
def dkl(s0,s1):
# modified Kullbaclk Lieber distance
#################################################################################### 
 m,n = s0.shape
 ds1 = LA.det(s1)
 ds0 = LA.det(s0)
 s1i = LA.inv(s1)
 a   = np.trace(np.matmul(s1i,s0)) - m 
 c   = np.log(ds1/ds0)
 dv  = 0.5*(a+c)
 return dv
  
####################################################################################   
def dinf(a,b):  
# infinite distance between a and b
#################################################################################### 
 m,n = a.shape
 vmax = 0
 for i in range(n):
  for j in range(i):
   vact = np.abs(a[i,j]-b[i,j]) 
   if vmax < vact:
    vmax = vact
 return vmax

####################################################################################     
def percentThres(U,pc):
####################################################################################   
 m,n = U.shape
 T = np.zeros((m,n))
 vx = triang(U)
 vpc = np.percentile(vx,pc)  
 for i in range(m):
  for j in range(n):
   if i == j:
    T[i,j] = 1
   else:
    if np.abs(U[i,j]) >= vpc:
     if U[i,j] < 0:
      T[i,j] = -1
     elif U[i,j] > 0:
      T[i,j] = 1
 return (T)


####################################################################################
def norml1(U):
####################################################################################
 t = triang(U)
 nl1 = np.sum(np.abs(t))
 return nl1

####################################################################### 
def diagones(A):
#######################################################################
 m,n =A.shape
 for i in range(m):
  for j in range(n):
   if i == j:
    A[i,j] = 1

####################################################################################     
def triang(U):
####################################################################################   
 m,n = U.shape
 vx = np.zeros(int((( n*n ) - n) / 2 ))
 k = 0
 for i in range(1,n):
  for j in range(i):
   vx[k] = U[i,j]
   k += 1
 return vx


####################################################################################      		 
def tpr(Q,X):
####################################################################################      		 
 tn = 0
 td = 0
 m,n  = Q.shape
 for i in range(m):
  for j in range(i):
   if X[i,j] == 0 and Q[i,j] == 0:
    tn += 1
   if Q[i,j] == 0:
    td += 1
 if td == 0:
  t = -1
 else:
  t = tn / td
 return t
 
####################################################################################      		 
def fpr(Q,X):
####################################################################################      		 
 fn = 0
 fd = 0
 m,n  = Q.shape
 for i in range(m):
  for j in range(i):
   if X[i,j] == 0 and Q[i,j] != 0:
    fn += 1
   if Q[i,j] != 0:
    fd += 1
 if fd == 0:
  f = -1
 else:
  f = fn / fd
 return f
 
  
 
####################################################################################      		 
def tnr(Q,X):
####################################################################################      		 
 tn = 0
 td = 0
 m,n  = Q.shape
 for i in range(m):
  for j in range(i):
   if X[i,j] != 0 and Q[i,j] != 0:
    tn += 1
   if Q[i,j] != 0:
    td += 1
 if td == 0:
  t = -1
 else:
  t = tn / td
 return t
 
####################################################################################      		 
def fnr(Q,X):
####################################################################################      		 
 fn = 0
 fd = 0
 m,n  = Q.shape
 for i in range(m):
  for j in range(i):
   if X[i,j] != 0 and Q[i,j] == 0:
    fn += 1
   if Q[i,j] == 0:
    fd += 1
 if fd == 0:
  f = -1
 else:
  f = fn / fd
 return f
 

###################################################################################
def tpro(Q,X):
###################################################################################
 tprn = 0
 tprd = 0
 m,n  = Q.shape
 for i in range(m):
  for j in range(i):
   if X[i,j] > 0 and Q[i,j] > 0:
    tprn += 1
   if Q[i,j] > 0:
    tprd += 1
 if tprd == 0:
  t = -1
 else:
  t = tprn / tprd
 return t
 
####################################################################################      		 
def fpro(Q,X):
####################################################################################      		 
 fprn = 0
 fprd = 0
 m,n  = Q.shape
 for i in range(m):
  for j in range(i):
   if X[i,j] > 0 and Q[i,j] == 0:
    fprn += 1
   if Q[i,j] == 0:
    fprd += 1
 if fprd == 0:
  f = -1
 else:
  f = fprn / fprd
 return f

####################################################################################
def defpath(vpath):
####################################################################################
 stop = len(vpath)-1
 k = len(vpath)
 for i in range(stop,0,-1):
  if vpath[i]=='/':
    k = i
    break
 path = vpath[:k]
 return path

