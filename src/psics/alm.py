# psics python package
# file: alm.py
# calls alm algorithm using call_alm.m
#
import numpy as np
import subprocess
import time
import importlib.resources as resources
from psics import utils

#######################################################################
def call_alm(*args, **kwargs):
#######################################################################

 if(len(args)<1):
   print("Error: ALM requires at least sigma matrix as input") 
   return

# 1) save covariance data to temp file
 sigma = args[0]
 np.savetxt('_almsigma.dat',sigma,delimiter=',',newline='\n')

# 2) process input parameters
 # default values
 rho = str(0.5)       # default rho value
 mxitr = str(500)     # max iteration number
 mu0 = str(1e-1)      # initial mu
 muf = str(1e-3)      # final mu
 rmu = str(1/4)       # ratio of decreasing mu
 tol_gap = str(1e-1)  # tolerance for duality gap
 tol_frel = str(1e-7) # tolerance for relative change of obj value
 tol_Xrel = str(1e-7) # tolerance for relative change of X
 tol_Yrel = str(1e-7) # tolerance for relative change of Y
 tol_pinf = str(1e-3) # tolerance for infeasibility
 numDG = str(10)      # every numDG iterations, we compute duality gap since it's expensive
 record = str(1)      # print stats
 fsigma = str(1e-10)  # fsigma is the smoothness parameter
 timeout = None
 for key,value in  kwargs.items():
   # rho
   if key == "rho":
    if (isinstance(value,(float))):
      rho = str(value)
	 # mxitr
   elif key== "mxitr":
    if (isinstance(value,(int,float))):
     mxitr = str(value)
	 # mu0
   elif key== "mu0":
    if (isinstance(value,(float))):
     mu0 = str(value)
	 # muf
   elif key== "muf":
    if (isinstance(value,(float))):
     muf = str(value)
	 # rmu
   elif key== "rmu":
    if (isinstance(value,(float))):
     rmu = str(value)
   # tol_gap
   elif key == "tol_gap":
    if(isinstance(value,(float))):
     tol_gap = str(value)
   # tol_frel
   elif key == "tol_frel":
    if(isinstance(value,(float))):
     tol_frel = str(value)
   # tol_Xrel
   elif key == "tol_Xrel":
    if(isinstance(value,(float))):
     tol_Xrel = str(value)
   # tol_Yrel
   elif key == "tol_Yrel":
    if(isinstance(value,(float))):
     tol_Yrel = str(value)
   # tol_pinf
   elif key == "tol_pinf":
    if(isinstance(value,(float))):
     tol_pinf = str(value)
   # numDG
   elif key == "numDG":
    if(isinstance(value,(float,int))):
     numDG = str(value)
   # record
   elif key == "record":
    if(isinstance(value,(float,int))):
     record = str(value)
   # fsigma
   elif key == "fsigma":
    if(isinstance(value,(float,int))):
     fsigma = str(value)
   elif key == "timeout":
    if(isinstance(value,(int,float))):
     timeout = value

 # 3) call system subprocess using octave

 # Get the file content and  path in a package
 ref = resources.files('psics') / 'call_alm.m'
 with resources.as_file(ref) as path:
   vpath = path
   path = utils.defpath(str(vpath))
 timestart = time.time()
 try:
   cp = subprocess.run(["octave", "--no-gui", str(vpath),"_almsigma.dat",
   rho,
   mxitr,
   mu0,
   muf,
   rmu,
   tol_gap,
   tol_frel,
   tol_Xrel,
   tol_Yrel,
   tol_pinf,
   numDG,
   record,
   fsigma,
   path
   ],
   timeout = timeout,
   check = True,
   capture_output = True,
   text = True)

   elaptime = time.time() - timestart

   # 4) read results files created by alm
   # inverse covariance, covariance estimated
   X = np.loadtxt('_almXsol.dat',dtype='float',delimiter=',')
   Y = np.loadtxt('_almYsol.dat',dtype='float',delimiter=',')

   # additional results
   r = np.loadtxt('_almres.dat',dtype='float',delimiter=',')

   niter = r[0]
   pinf  = r[1]
   obj   = r[2]
   gapX  = r[3]
   gapY  = r[4]
   gap   = r[5]

   # 5) delete support files
   df = subprocess.run("rm _*.dat", shell=True)

   return {
     'X':            X,
     'Y':            Y,
     'niter':        niter,
     'pinf':         pinf,
     'obj':          obj,
     'gapX':         gapX,
     'gapY':         gapY,
     'gap':          gap,
     'success' :     True,
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

