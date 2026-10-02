# psics python package
# file: covsel.py
# calls covsel algorithm using call_covsel.m
#
import numpy as np
import subprocess
import time
import importlib_resources as resources
from psics import utils

#######################################################################
def call_covsel(*args, **kwargs):
#######################################################################
#
# 1 save covariance matrix in auxiliar data file

 if(len(args)<1):
  print("Error: COVSEL requires at least sigma matrix as input")
  return

# 1) save covariance data to temp file
 sigma = args[0]
 np.savetxt('_covselsigma.dat',sigma,delimiter=',',newline='\n')

# 2) process input parameters
 # default values
 rho     = str(0.5)       # default rho value
 maxiter = str(100)     # max iteration number
 prec    = str(1e-1)
 maxnest = str(1e2)
 algo    = "BoxQP"
 timeout = None

 for key,value in  kwargs.items():
  # rho
  if key == "rho":
   if (isinstance(value,(float))):
    rho = str(value)
	# mxitr
  elif key== "maxiter":
   if (isinstance(value,(int,float))):
    maxiter = str(value)
	# mu0
  elif key== "prec":
   if (isinstance(value,(float))):
    prec = str(value)
	# muf
  elif key == "maxnest":
   if (isinstance(value,(int,float))):
    maxnest = str(value)
	# algot
  elif key == "algo":
   if (key == "BoxQP" or key == "sedumi"):
    algot = value

  elif key == "timeout":
   if(isinstance(value,(int,float))):
    timeout = value

# 3) call system subprocess using octave

 ref = resources.files('psics') / 'call_covsel.m'
 with resources.as_file(ref) as path:
  vpath = path
 path = utils.defpath(str(vpath))
 timestart = time.time()
 try:
   cp = subprocess.run(["octave", "--no-gui",str(vpath),"_covselsigma.dat",
        rho,
        maxiter,
        prec,
        maxnest,
        algo,
        path],
        timeout=timeout,
        check=True,
        text=True,
        capture_output=True)

   elaptime = time.time() - timestart

# 4) read results
# inverse covariance, covariance estimated and relative gap
   X = np.loadtxt('_covselXsol.dat',dtype='float',delimiter=',')
   U = np.loadtxt('_covselUsol.dat',dtype='float',delimiter=',')
   gvals =  np.loadtxt('_covselgvals.dat',dtype='float',delimiter=',')
   cputimes = np.loadtxt('_covselcputimes.dat',dtype='float',delimiter=',')

# 5) delete support files
   df = subprocess.run("rm _*.dat", shell=True)

   return {
     'X':            X,
     'U':            U,
     'gvals':        gvals,
     'cputimes':     cputimes,
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
     'stdout' :      None,
     'stderr' :      None,
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


