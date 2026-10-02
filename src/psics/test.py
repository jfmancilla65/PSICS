# psics python package
# file: test.py
# python function to test sics algorithms covsel, covpath, alm, glasso and glassopath
# creates a synthetic covariance matrix as input
# generates figures and table of numerical results
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from   psics import covpath,utils,covsel,glasso,alm
import numpy.linalg as LA
import math as ma
import os

#----------------------------------------------------------------------
# examples using  covariance selection algorithms
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def testall(inst):
#----------------------------------------------------------------------
# function for testing all algorithms for the instance defined in inst
# if inst == "ld", low dimension instance is tested
# if inst == "hd", high dimension instance is tested

 # synthetic test data
 eps = 1e-4
 # dimension matrix
 n = 60
 pc = 92
 # number of Gaussian samples
 if inst == "ld" or inst.upper()=="LD":
  ns = 5 * n      # low dimension
 else:
  ns = int(n / 4) # high dimension

 rng = np.random.default_rng(n)

 # ratio of nonzero elements of sparse matrix
 nnz = 10 / n
 ne = int(((n * n)-n) / 2 * nnz )   # elements not in main diagonal
 A = np.eye(n)

 # synthetic sparse inverse covariance matrix
 for i in range(ne):
  im = 0
  jm = 0
  while im == jm or A[im,jm] != 0:
   im = ma.floor(n*rng.random())
   jm = ma.floor(n*rng.random())
  nm = np.sign(rng.random()-.5)
  A[im,jm] = nm
  A[jm,im] = nm
 d = np.diag(A).copy()
 D = np.diag(d)
 for i in range(n):
  for j in range(n):
   A[i,j] = min(A[i,j]-D[i,j],1)
   A[i,j] = max(A[i,j],-1)
 A = A + np.diag(d+1)

 # make matrix positive definite
 U = A + max(-1.2*min(LA.eigvalsh(A)),eps)*np.eye(n)

 # B is the ground-truth covariance matrix
 B = LA.inv(U)
 data = rng.multivariate_normal(np.zeros(n),B,ns)
 S = (1/n)*np.matmul(np.transpose(data),data)
 S = S + max(-1.2*min(LA.eigvalsh(S)),eps)*np.eye(n)


#----------------------------------------------------------------------
# covpath
#----------------------------------------------------------------------
 print("COVPATH test")
 m,n = S.shape

 # rholist vector
 rhomax = np.max(np.diag(S))
 nrhos = 10
 logrhomax = np.log10(rhomax)
 logrhomin = logrhomax - 1.5
 rholist = 10 ** (np.linspace(logrhomax, logrhomin, nrhos + 1))
 # sort rholist (decreasing)
 rholist = np.sort(rholist)[::-1]
 rhomax = 1.2 * rhomax
 outcovpath  = covpath.call_covpath(S,rholist=rholist,t=1e-5,colFraction=1,rhomax=rhomax)
 X   = outcovpath.get('invCovmat') # inverse sparse matrices for all values of rholist
 Y   = outcovpath.get('Covmat')    # covariance calculated matrices for all values of rholist
 Xcovpath = covpath.covpath_matrix(X,m,9) # read inverse sparse matrix for rholist[9]
 Ycovpath = covpath.covpath_matrix(Y,m,9) # read covariance matrix for rholist[9]
 Ycovpath = Ycovpath + max(-1.2*min(LA.eigvals(Ycovpath)),eps)*np.eye(n)
 print("COVPATH output")
 print("Xcovpath=")
 print(Xcovpath)
 print("Ycovpath=")
 print(Ycovpath)
 print('rholist= ',outcovpath.get('rholist'))
 print('relative gaps= ',outcovpath.get('relativeGaps'))
 print("returncode= ", outcovpath.get('returncode'))
 print("shell exec time= ", outcovpath.get('elaptime'))
 Xcovpath_pc = utils.percentThres(Xcovpath,pc)
 utils.diagones(Xcovpath_pc)
 resultstab = [("COVPATH",utils.d1(U,Xcovpath_pc),utils.d2(U,Xcovpath_pc),
                         utils.dinf(U,Xcovpath_pc),utils.dkl(S,Ycovpath),
                         utils.norml1(Xcovpath_pc),
                         utils.tpr(U,Xcovpath_pc),utils.fpr(U,Xcovpath_pc),
                         utils.tnr(U,Xcovpath_pc),utils.fnr(U,Xcovpath_pc),
                         utils.tpro(U,Xcovpath_pc),utils.fpro(U,Xcovpath_pc))]

#----------------------------------------------------------------------
# covsel
#----------------------------------------------------------------------
 print("COVSEL test")
 # call COVSEL algorithm
 rhocovsel = 0.01
 outcovsel = covsel.call_covsel(S,rho=rhocovsel, maxiter=100,prec=1.e-1,maxnest=100,algot='BoxQP')
 print("COVSEL output")
 Xcovsel = outcovsel.get('X')
 Ycovsel = LA.inv(Xcovsel)
 print("Xcovsel=")
 print(Xcovsel)
 print("Ycovsel=")
 print(Ycovsel)
 print("gvals= ",outcovsel.get('gvals'))
 print("cputimes= ",outcovsel.get('cputimes'))
 print("returncode= ",outcovsel.get('returncode'))
 print("shell exec time= ", outcovsel.get('elaptime'))
 Ycovsel = Ycovsel + max(-1.2*min(LA.eigvals(Ycovsel)),eps)*np.eye(n)
 Xcovsel_pc = utils.percentThres(Xcovsel,pc)
 utils.diagones(Xcovsel_pc)
 resultstab.append(("COVSEL",utils.d1(U,Xcovsel_pc),utils.d2(U,Xcovsel_pc),
                            utils.dinf(U,Xcovsel_pc),utils.dkl(S,Ycovsel),
                            utils.norml1(Xcovsel_pc),
                            utils.tpr(U,Xcovsel_pc),utils.fpr(U,Xcovsel_pc),
                            utils.tnr(U,Xcovsel_pc),utils.fnr(U,Xcovsel_pc),
                            utils.tpro(U,Xcovsel_pc),utils.fpro(U,Xcovsel_pc)))


#----------------------------------------------------------------------
# alm
#----------------------------------------------------------------------
 print("ALM test")
 # call ALM algorithm
 rhoalm = 0.01
 mu0 = 100/rhoalm
 outalm = alm.call_alm(S,mxitr=100,mu0=mu0,muf=0.0001,rmu=0.25,tol_gap=1e-3,tol_frel=1e-7,
                tol_Xrel=1e-7,tol_Yrel=1e-7,tol_pinf=1e-1,numDG=10,record=1,fsigma=1e2)
 print("ALM output")
 Xalm = outalm.get("X")
 print("Xalm=")
 print(Xalm)
 Yalm = LA.inv(Xalm)
 Yalm = Yalm + max(-1.2*min(LA.eigvals(Yalm)),eps)*np.eye(n)
 print("Yalm=")
 print(Yalm)
 print("pinf= ",outalm.get("pinf"))
 print("obj = ",outalm.get("obj"))
 print("gapX= ",outalm.get("gapX"))
 print("gapY= ",outalm.get("gapY"))
 print("gap = ",outalm.get("gap"))
 print("returncode= ", outalm.get('returncode'))
 print("shell exec time= ", outalm.get('elaptime'))

 Xalm_pc = utils.percentThres(Xalm,pc)
 utils.diagones(Xalm_pc)
 resultstab.append(("ALM",utils.d1(U,Xalm_pc),utils.d2(U,Xalm_pc),
                         utils.dinf(U,Xalm_pc),utils.dkl(S,Yalm),
                         utils.norml1(Xalm_pc),
                         utils.tpr(U,Xalm_pc),utils.fpr(U,Xalm_pc),
                         utils.tnr(U,Xalm_pc),utils.fnr(U,Xalm_pc),
                         utils.tpro(U,Xalm_pc),utils.fpro(U,Xalm_pc)))

#----------------------------------------------------------------------
# glasso
#----------------------------------------------------------------------
 print("GLASSO test")
 # call GLASSO algorithm
 rho = 0.01
 outglasso = glasso.call_glasso(S,rho=rho,trace=True,approx=False,maxit=1e5,thr=1e-4)
 Xglasso = outglasso.get("wi")
 Yglasso = outglasso.get("w")
 errflag = outglasso.get("errflag")
 print("GLASSO output")
 print("Xglasso=")
 print(Xglasso)
 print("Yglasso=")
 print(Yglasso)
 print("errflag = ",errflag)
 print("returncode= ", outglasso.get('returncode'))
 print("shell exec time=  ", outglasso.get('elaptime'))
 Yglasso = Yglasso + max(-1.2*min(LA.eigvals(Yglasso)),eps)*np.eye(n)
 Xglasso_pc = utils.percentThres(Xglasso,pc)
 utils.diagones(Xglasso_pc)
 resultstab.append(("GLASSO",utils.d1(U,Xglasso_pc),utils.d2(U,Xglasso_pc),
                            utils.dinf(U,Xglasso_pc),utils.dkl(S,Yglasso),
                            utils.norml1(Xglasso_pc),
                            utils.tpr(U,Xglasso_pc),utils.fpr(U,Xglasso_pc),
                            utils.tnr(U,Xglasso_pc),utils.fnr(U,Xglasso_pc),
                            utils.tpro(U,Xglasso_pc),utils.fpro(U,Xglasso_pc)))

#----------------------------------------------------------------------
# glassopath
#----------------------------------------------------------------------
 print("GLASSOPATH test")
 # call GLASSOPATH algorithm
 m,n = S.shape

 # rholist vector
 rhomax = np.max(np.diag(S))
 nrhos = 10
 logrhomax = np.log10(rhomax)
 logrhomin = logrhomax - 1.5
 rholist = 10 ** (np.linspace(logrhomax, logrhomin, nrhos + 1))
 # sort rholist
 rholist = np.sort(rholist)

 outglassopath = glasso.call_glassopath(S,rholist=rholist,trace=1,thr=1e-5)
 X = outglassopath.get("wi") # inverse sparse matrices for all values of rholist
 Y = outglassopath.get("w")  # covariance calculated matrices for all values of rholist
 Xglassopath = glasso.glassopath_matrix(X,m,1) # read inverse sparse matrix for rholist[1]
 Yglassopath = glasso.glassopath_matrix(Y,m,1) # read covariance matrix for rholist[1]
 Yglassopath = Yglassopath + max(-1.2*min(LA.eigvals(Yglassopath)),eps)*np.eye(n)
 print("GLASSOPATH output")
 print("Xglassopath=")
 print(Xglassopath)
 print("Yglassopath=")
 print(Yglassopath)
 errflag = outglassopath.get("errflag")
 rl =      outglassopath.get("rholist")
 print("errflag path= ",errflag)
 print("rholist= ",rl)
 print("returncode= ", outglassopath.get('returncode'))
 print("shell exec time= ", outglassopath.get('elaptime'))
 Xglassopath_pc = utils.percentThres(Xglassopath,pc)
 utils.diagones(Xglassopath_pc)
 resultstab.append(("GLASSOPATH", utils.d1(U,Xglassopath_pc),utils.d2(U,Xglassopath_pc),
                                 utils.dinf(U,Xglassopath_pc),utils.dkl(S,Yglassopath),
                                 utils.norml1(Xglassopath_pc),
                                 utils.tpr(U,Xglassopath_pc),utils.fpr(U,Xglassopath_pc),
                                 utils.tnr(U,Xglassopath_pc),utils.fnr(U,Xglassopath_pc),
                                 utils.tpro(U,Xglassopath_pc),utils.fpro(U,Xglassopath_pc)))

#----------------------------------------------------------------------
# print results
#----------------------------------------------------------------------

 # distance/divergence measures
 dklreal = utils.dkl(S,S)

 # additional dkl tests
 ide = np.identity(m)
 dia = np.diag(np.diag(U))
 dri = utils.dkl(U,ide)
 drd = utils.dkl(U,dia)
 dri = utils.dkl(ide,U)
 drd = utils.dkl(dia,U,)
 print("Instance test = ",inst)
 print("dkl(U,ide) = ",dri)
 print("dkl(U,diag) =  ",drd)
 print("dkl(ide,U) = ",dri)
 print("dkl(diag,U) =  ",drd)
 print("norml1(U) = ", utils.norml1(U))

 print("Algorithm\t DL1\t\t DL2\t\t DINF\t\t DKL\t\t |X|1\t\tTPR\tFPR\tTNR\tFNR\tTPRO\tFPRO") 
 for i, (name, d1, d2, dinf, dkl, nl1,tpr,fpr,tnr,fnr,tpro,fpro) in enumerate(resultstab):
  print(f"{name:15s}\t{d1:8.3f}\t{d2:8.3f}\t{dinf:8.3f}\t{dkl:8.3f}\t{nl1:8.3f}\t{tpr:.3f}\t{fpr:.3f}\t{tnr:.3f}\t{fnr:.3f}\t{tpro:.3f}\t{fpro:.3f}")


#----------------------------------------------------------------------
# plot inverses
#----------------------------------------------------------------------
 utils.diagones(U)
 cmap = mpl.colors.ListedColormap(['orange', 'white', 'blue'])

 # find minimum of minima & maximum of maxima
 minmin = np.min([np.min(U), np.min(Xcovpath_pc),np.min(Xalm_pc),np.min(Xcovsel_pc),np.min(Xglasso_pc),np.min(Xglassopath_pc)])
 maxmax = np.max([np.max(U), np.max(Xcovpath_pc),np.max(Xalm_pc),np.max(Xcovsel_pc),np.max(Xglasso_pc),np.max(Xglassopath_pc)])

 # all matrices
 fig, axes = plt.subplots(nrows=2, ncols=3)
 axes[0][0].set_title("Real inverse",fontsize='small',family='monospace')
 axes[0][0].tick_params(axis='both',labelsize='small')
 axes[0][1].set_title("COVSEL",fontsize='small',family='monospace')
 axes[0][1].tick_params(axis='both',labelsize='small')
 axes[0][2].set_title("ALM",fontsize='small',family='monospace')
 axes[0][2].tick_params(axis='both',labelsize='small')
 axes[1][0].set_title("GLASSO",fontsize='small',family='monospace')
 axes[1][0].tick_params(axis='both',labelsize='small')
 axes[1][1].set_title("COVPATH",fontsize='small',family='monospace')
 axes[1][1].tick_params(axis='both',labelsize='small')
 axes[1][2].set_title("GLASSOPATH",fontsize='small',family='monospace')
 axes[1][2].tick_params(axis='both',labelsize='small')

 im1 = axes[0,0].imshow(U, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im2 = axes[0,1].imshow(Xcovsel_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im3 = axes[0,2].imshow(Xalm_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im4 = axes[1,0].imshow(Xglasso_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im5 = axes[1,1].imshow(Xcovpath_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im6 = axes[1,2].imshow(Xglassopath_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 fig.subplots_adjust(right=0.85)
 cbar_ax = fig.add_axes([0.9, 0.27, 0.005, 0.45]) 
 cbar_ax.tick_params(labelsize='small')
 fig.colorbar(im2, cax=cbar_ax,ticks=[minmin,0,maxmax])
 plt.savefig("sicsmatrix_all_"+inst+".jpg",dpi=800)
 plt.close()



 # covpath
 fig, axes = plt.subplots(nrows=1, ncols=2)
 axes[0].set_title("Real inverse",fontsize='small',family='monospace')
 axes[1].set_title("COVPATH",fontsize='small',family='monospace')
 axes[0].tick_params(axis='both',labelsize='small')
 axes[1].tick_params(axis='both',labelsize='small')
 im1 = axes[0].imshow(U, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im2 = axes[1].imshow(Xcovpath_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 fig.subplots_adjust(right=0.85)
 cbar_ax = fig.add_axes([0.9, 0.27, 0.005, 0.45]) 
 cbar_ax.tick_params(labelsize='small')
 fig.colorbar(im2, cax=cbar_ax,ticks=[minmin,0,maxmax])
 plt.savefig("sicsmatrix_COVPATH_"+inst+".jpg",dpi=800)
 plt.close()

 # covsel
 fig, axes = plt.subplots(nrows=1, ncols=2)
 axes[0].set_title("Real inverse",fontsize='small',family='monospace')
 axes[1].set_title("COVSEL",fontsize='small',family='monospace')
 axes[0].tick_params(axis='both',labelsize='small')
 axes[1].tick_params(axis='both',labelsize='small')
 im1 = axes[0].imshow(U, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im2 = axes[1].imshow(Xcovsel_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 fig.subplots_adjust(right=0.85)
 cbar_ax = fig.add_axes([0.9, 0.27, 0.005, 0.45]) 
 cbar_ax.tick_params(labelsize='small')
 fig.colorbar(im2, cax=cbar_ax,ticks=[-1,0,1])
 plt.savefig("sicsmatrix_COVSEL_"+inst+".jpg",dpi=800)
 plt.close()

 # alm
 fig, axes = plt.subplots(nrows=1, ncols=2)
 axes[0].set_title("Real inverse",fontsize='small',family='monospace')
 axes[1].set_title("ALM",fontsize='small',family='monospace')
 axes[0].tick_params(axis='both',labelsize='small')
 axes[1].tick_params(axis='both',labelsize='small')
 im1 = axes[0].imshow(U, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im2 = axes[1].imshow(Xalm_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 fig.subplots_adjust(right=0.85)
 cbar_ax = fig.add_axes([0.9, 0.27, 0.005, 0.45]) 
 cbar_ax.tick_params(labelsize='small')
 fig.colorbar(im2, cax=cbar_ax,ticks=[-1,0,1])
 plt.savefig("sicsmatrix_ALM_"+inst+".jpg",dpi=800)
 plt.close()

 # glasso
 fig, axes = plt.subplots(nrows=1, ncols=2)
 axes[0].set_title("Real inverse",fontsize='small',family='monospace')
 axes[1].set_title("GLASSO",fontsize='small',family='monospace')
 axes[0].tick_params(axis='both',labelsize='small')
 axes[1].tick_params(axis='both',labelsize='small')
 im1 = axes[0].imshow(U, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im2 = axes[1].imshow(Xglasso_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 fig.subplots_adjust(right=0.85)
 cbar_ax = fig.add_axes([0.9, 0.27, 0.005, 0.45]) 
 cbar_ax.tick_params(labelsize='small')
 fig.colorbar(im2, cax=cbar_ax,ticks=[-1,0,1])
 plt.savefig("sicsmatrix_GLASSO_"+inst+".jpg",dpi=800)
 plt.close()

 # glassopath
 fig, axes = plt.subplots(nrows=1, ncols=2)
 axes[0].set_title("Real inverse",fontsize='small',family='monospace')
 axes[1].set_title("GLASSOPATH",fontsize='small',family='monospace')
 axes[0].tick_params(axis='both',labelsize='small')
 axes[1].tick_params(axis='both',labelsize='small')
 im1 = axes[0].imshow(U, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 im2 = axes[1].imshow(Xglassopath_pc, vmin=minmin, vmax=maxmax,cmap=cmap,extent=(0,n,0,n))
 fig.subplots_adjust(right=0.85)
 cbar_ax = fig.add_axes([0.9, 0.27, 0.005, 0.45]) 
 cbar_ax.tick_params(labelsize='small')
 fig.colorbar(im2, cax=cbar_ax,ticks=[-1,0,1])
 plt.savefig("sicsmatrix_GLASSOPATH_"+inst+".jpg",dpi=800)
 plt.close()


#--------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------

def run():
 # set R environment variable to access the user library for psics
 psicsRlib = os.environ.get("HOME")
 psicsRlib = psicsRlib + "/R/lib/psics" 
 if not os.environ.get("R_LIBS_USER"):
  os.environ.setdefault("R_LIBS_USER",psicsRlib)
 else:
  os.environ["R_LIBS_USER"] = psicsRlib

 # ----------------------------------------------
 # run sics algorithms demo
 # ----------------------------------------------

 # low dimension instance
 testall("ld")

 # high dimension instance
 testall("hd")



