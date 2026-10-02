# psics python package
# file: glasso.py
# calls glasso and glasso path algorithms using call_glasso.R
# and call_glassopathR.R
#
import numpy as np
import subprocess
import time
import importlib.resources

#######################################################################
def call_glassopath(*args, **kwargs):
#######################################################################

 if(len(args)<1):
  print("Error: Glassopath requires at least sigma matrix as input") 
  return 

# arguments default values

 rholist = "NULL"
 thr     = str(1.0e-4)
 maxit   = str(1e4)
 approx  = "FALSE"
 penalizediagonal = "TRUE"
 winit  = "NULL"
 wiinit = "NULL"
 trace   = str(0)
 timeout = None

# 1) save covariance data to temp file
 sigma = args[0]
 np.savetxt('_glassosigma.dat',sigma,delimiter=',',newline='\n')

# 2) process input parameters
# default parameters values

 for key,value in  kwargs.items():
  if key == "rholist":
   np.savetxt("_glassorho.dat",value,delimiter=',')
   rholist = "_glassorho.dat"
	# thr
  elif key== "thr":
   if (isinstance(value,(float))):
    thr = str(value)
	# maxit
  elif key== "maxit":
   if (isinstance(value,(int,float))):
    maxit = str(value)
	# approx
  elif key== "approx":
   if value:
    approx = "TRUE"
	# penalizediagonal
  elif key== "penalizediagonal":
   if not value:
    penalizediagonal="FALSE"
     # winit
  elif key == "winit":
   winit = "_glassowinit.dat"
   np.savetxt(winit,value,delimiter=',')
    # wiinit
  elif key == "wiinit":
   wiinit = "_glassowiinit.dat"
   np.savetxt(wiinit,value,delimiter=',')
    # trace
  elif key == "trace":
   if(value >= 0 and value <= 2 ):
    trace = str(value)
  elif key == "timeout":
   if(isinstance(value,(int,float))):
    timeout = value
 # 5) call glassopath function
 ref = importlib.resources.files('psics') / 'call_glassopath.R'
 with importlib.resources.as_file(ref) as path:
  vpath = path
 timestart = time.time()
 try:
   cp = subprocess.run(["Rscript",str(vpath),"_glassosigma.dat",
   rholist,
   thr,
   maxit,
   approx,
   penalizediagonal,
   winit,
   wiinit,
   trace
   ],
   timeout=timeout,
   capture_output=True,
   check=True,
   text=True)

   elaptime = time.time()-timestart

# 3) read results from temp files
#   glasso rholist 
   rl = np.loadtxt('_glassopathrholist.dat',dtype='float',delimiter=',')
#   error list
   el = np.loadtxt('_glassopatherrflag.dat',dtype='float',delimiter=',')
#   inverse covariance and covariance estimated
   X = np.loadtxt('_glassopathXsol.dat',dtype='float',delimiter=',')
   U = np.loadtxt('_glassopathUsol.dat',dtype='float',delimiter=',')
# other results defined by glasso
   with open('_glassopathres.dat') as file:
    resources = {}
    for line in file:
     key, value = line.rstrip().split(',', 1)
     resources[key] = value
   approx  = resources['approx']

# 4) delete support files

   df = subprocess.run("rm _*.dat", shell=True)

   return {
     'w':   U,
     'wi':X,
     'approx':       approx,
     'rholist':      rl,
     'errflag':      el,
     'elaptime':     elaptime,
     'success':      True,
     'returncode':   cp.returncode,
     'stdout' :      cp.stdout.strip(),
     'stderr' :      cp.stderr.strip(),
     'errortype' :   None,
     'elaptime' :    elaptime
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
def glassopath_matrix(v,n,ind):
####################################################################################
# Selects ind matrix from v matrix covariance or covariance inverse results.
# All path partial result matrices were stored using rholist value as row
# Each row contains a n*n matriz stored by rows
# 
 mat = np.zeros((n,n))
 m,p = v.shape
 if ind > m or ind < 0:
  return
 k = 0
 for i in range(n):
  for j in range(n):
   mat[i,j] = v[ind-1,k] 
   k += 1
 return mat

 
#######################################################################
def call_glasso(*args, **kwargs):
#######################################################################
 
 if(len(args)<1):
  print("Error: Glasso requires at least sigma matrix as input") 
  return 

# 1) save covariance data to temp file
 sigma = args[0]
 np.savetxt('_glassosigma.dat',sigma,delimiter=',',newline='\n')

# 2) arguments default values

 nobs    = "NULL"
 zero    = "NULL"
 thr     = str(1.0e-4)
 maxit   = str(1e4)
 approx  = "FALSE"
 penalizediagonal = "TRUE"
 start   = "cold"
 winit   = "NULL"
 wiinit  = "NULL"
 zero    = "NULL"
 trace   = "FALSE"
 rho     = str(0.5)
 timeout = None

# 3) process input parameters
 for key,value in  kwargs.items():
  # rho
  if key == "rho":
   if (isinstance(value,(int,float))):
    rho = str(value)
   else:
    np.savetxt("_glassorho.dat",rho,delimiter=',',newline='\n')
    rho = "_glassorho.dat"
  # nobs
  elif key == "nobs":
   if (isinstance(value,(int,float))):
    nobs = str(value)
  # zero
  elif key == "zero":
   zero = "_glassozeros.dat"
   np.savetxt(zero,value,delimiter=',')
  # thr
  elif key== "thr":
   if (isinstance(value,(float))):
    thr = str(value)
  # maxit
  elif key== "maxit":
   if (isinstance(value,(int,float))):
    maxit = str(value)
  # approx
  elif key== "approx":
   if value:
    approx = "TRUE"
  # penalizediagonal
  elif key== "penalizediagonal":
   if not value:
    penalizediagonal = "FALSE"
  # start
  elif key == "start":
   if(value == "warm" or value == "WARM"):
    start = "warm"
  # winit
  elif key == "winit":
   winit = "_glassowinit.dat"
   np.savetxt(winit,value,delimiter=',')
  # wiinit
  elif key == "wiinit":
   wiinit = "_glassowiinit.dat"
   np.savetxt(wiinit,value,delimiter=',')
  # trace
  elif key == "trace":
   if(value):
    trace = "TRUE"
  elif key == "timeout":
   if(isinstance(value,(int,float))):
    timeout = value

 # 4) call glasso function
 ref = importlib.resources.files("psics") / "call_glasso.R"
 with importlib.resources.as_file(ref) as path:
  vpath = path
 timestart = time.time()
 try:
   cp = subprocess.run(["Rscript", str(vpath),"_glassosigma.dat",
        rho,
        nobs,
        zero,
        thr,
        maxit,
        approx,
        penalizediagonal,
        start,
        winit,
        wiinit,
        zero,
        trace
        ],
        timeout=timeout,
        capture_output=True,
        check=True,
        text=True)

   elaptime = time.time() - timestart

# 5) read results from temp files
#    inverse covariance and covariance estimated
   X = np.loadtxt('_glassoXsol.dat',dtype='float',delimiter=',')
   U = np.loadtxt('_glassoUsol.dat',dtype='float',delimiter=',')

# other results defined by glasso
   with open('_glassores.dat') as file:
    resources = {}
    for line in file:
     key, value = line.rstrip().split(',', 1)
     resources[key] = value
   loglik  = resources['loglik']
   errflag = resources['errflag']
   approx  = resources['approx']
   vdel    = resources['del']
   niter   = resources['niter']

# 6) delete support files
   df = subprocess.run("rm _*.dat", shell=True)
   return {
     'wi':           X,
     'w':            U,
     'loglik':       loglik,
     'errflag':      errflag,
     'approx':       approx,
     'del':          vdel,
     'niter':        niter,
     'returncode':   cp.returncode,
     'success':      True,
     'stdout' :      cp.stdout.strip(),
     'stderr' :      cp.stderr.strip(),
     'errortype' :   None,
     'elaptime' :    elaptime
   }
 except subprocess.CalledProcessError as e:
   return {
     'success' :     False,
     'returncode' :  e.returncode,
     'stdout' :      e.stdout.strip() if e.stdout else "",
     'stderr' :      e.stderr.strip() if e.stdout else "", 
     'errortype' :   'CalledProcessError',
     'elaptime' :    time.time() - timestart
   }

 except FileNotFoundError:
   return {
     'success' :     False,
     'returncode' :  None,
     'stdout' :      "",
     'stderr' :      "",
     'errortype' :   'FileNotFoundError',
     'elaptime' :    time.time() - timestart
   }

 except subprocess.TimeoutExpired as e:
   return {
     'success' :     False,
     'returncode' :  None,
     'stdout' :      e.stdout.strip() if e.stdout else "",
     "stderr" :      e.stderr.strip() if e.stdout else "",
     "errortype" :   "TimeoutExpired",
     "elaptime" :    time.time() - timestart
   }

 except Exception as e:
   return {
     "success" :     False,
     "returncode" :  None,
     "stdout" :      e.stdout.strip() if e.stdout else "",
     "stderr" :      str(e),
     "errortype" :   "UnexpectedError",
     "elaptime" :    time.time() - timestart
   }
