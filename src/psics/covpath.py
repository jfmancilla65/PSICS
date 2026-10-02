# psics python package
# file: covpath.py
# calls covpath algorithm using call_covpath.R
#
import numpy as np
import subprocess
import time
import importlib.resources as resources

#######################################################################
def call_covpath(*args,**kwargs):
#######################################################################
 if(len(args)<1):
  print("Error: COVPATH requires at least sigma matrix as input")
  return

# 1) save covariance data to temp file
 sigma = args[0]
 np.savetxt('_covpathsigma.dat',sigma,delimiter=',',newline='\n')

# 2) process input parameters
 # default values
 rholist = "NULL"
 colFraction = str(1)
 t = str(10^-7)
 timeout = None
 rhomax = "NULL"

 for key,value in kwargs.items():
  # rholist
  if key == "rholist":
   rholist = '_covpathrholist.dat'
   np.savetxt(rholist,value,delimiter=',',newline='\n')
  # rhomax
  elif key == "rhomax":
   if (isinstance(value,(int,float))):
    rhomax = str(value)
 	# t
  elif key== "t":
   if (isinstance(value,(int,float))):
    t = str(value)
	# colFraction
  elif key== "colFraction":
   if (isinstance(value,(int,float))):
    colFraction = str(value)
  # timeout
  elif key == "timeout":
   if(isinstance(value,(int,float))):
    timeout = value

# 3) call covpath function
 ref = resources.files("psics") / "call_covpath.R"
 with resources.as_file(ref) as path:
  vpath = path
 timestart = time.time()
 try:
   cp = subprocess.run(["Rscript",str(vpath),"_covpathsigma.dat",
   rholist,
   rhomax,
   t,
   colFraction],
   timeout=timeout,
   capture_output=True,
   check=True,
   text=True)

   elaptime = time.time() - timestart

# 4) read results files created by covpath
# inverse covariance, covariance estimated and relative gap
   X = np.loadtxt('_covpathXsol.dat',dtype='float',delimiter=',')
   U = np.loadtxt('_covpathUsol.dat',dtype='float',delimiter=',')
   rgap = np.loadtxt('_covpathRelGap.dat',dtype='float',delimiter=',')
   rholist = np.loadtxt('_covpathRholist.dat',dtype='float',delimiter=',')

# 5 delete support files
   df = subprocess.run("rm _*.dat", shell=True)

   return {'invCovmat':    X,
           'Covmat':       U,
           'relativeGaps': rgap,
           'rholist':      rholist,
           'elaptime':     elaptime,
           'returncode':   cp.returncode,
           'success':      True,
           'returncode':   cp.returncode,
           'stdout':       cp.stdout,
           'stderr':       cp.stderr,
           'errortype':    None,
           'elaptime':     elaptime
   }
 except subprocess.CalledProcessError as e:
   return {
     'success' :     False,
     'returncode' :  e.returncode,
     'stdout' :      e.stdout.strip() if e.stdout else "",
     'stderr' :      e.stderr.strip() if e.stdout else "", 
     'errortype' :  'CalledProcessError',
     'elaptime' : time.time() - timestart
   }

 except FileNotFoundError:
   return {
     'success' :     False,
     'returncode' :  None,
     'stdout' :      "",
     'stderr' :      "",
     'errortype' :  'FileNotFoundError',
     'elaptime' : time.time() - timestart
   }

 except subprocess.TimeoutExpired as e:
   return {
     'success' :     False,
     'returncode' :  None,
     'stdout' :      e.stdout.strip() if e.stdout else "",
     "stderr" :      e.stderr.strip() if e.stdout else "",
     "errortype" :  "TimeoutExpired",
     "elaptime" : time.time() - timestart
   }

 except Exception as e:
   return {
     "success" :     False,
     "returncode" :  None,
     "stdout" :      e.stdout.strip() if e.stdout else "",
     "stderr" :      str(e),
     "errortype" :  "UnexpectedError",
     "elaptime" : time.time() - timestart
   }

####################################################################################
def covpath_matrix(v,n,ind):
####################################################################################
# Selects ind matrix from v covpath vector covariance or covariance inverse results.
# All path partial result matrices were stored by row in vertical vector v.
#
 mat = np.zeros((n,n))

 k = n * n * (ind-1)
 for i in range(n):
  for j in range(n):
   mat[i,j] = v[k]
   k += 1
 return mat


