# psics python package
# file: config.py
# configure sics algorthms in python package
#
import zipfile
import tarfile
from urllib.request import urlretrieve
import re
import subprocess
import importlib_resources as resources
import os

from psics import utils

#----------------------------------------------------------
def modify(filepath, from_, to_):
#----------------------------------------------------------
 file = open(filepath,"r+")
 text = file.read()
 pattern = from_
 splitted_text = re.split(pattern,text)
 modified_text = to_.join(splitted_text)
 with open(filepath, 'w') as file:
  file.write(modified_text)

#----------------------------------------------------------
def setup():
#----------------------------------------------------------
 home = os.environ['HOME'] 
 dwnCovsel  = False
 urlCovsel  = "https://www.di.ens.fr/~aspremon/ZIP/COVSEL.zip"
 dwnCovpath = False
 urlCovpath = "https://cran.r-project.org/src/contrib/Archive/Covpath/Covpath_1.0.tar.gz"
 dwnAlm     = False
 urlAlm     = "https://github.com/linboqiao/SICS/archive/refs/heads/master.zip"
 urlGlasso  = "https://cran.r-project.org"

# define R library packages  in the user home directory
 psicsRlib  = home + "/R/lib/psics"

# ----------------------------------------------
# setup sics algorithms code
# ----------------------------------------------
 setupcode(dwnCovsel,dwnCovpath,dwnAlm,urlCovsel,urlCovpath,urlAlm,urlGlasso,psicsRlib)

 
#-------------------------------------------------------------------------------------------
def setupcode(dwnCovsel,dwnCovpath,dwnAlm,urlCovsel,urlCovpath,urlAlm,urlGlasso,psicsRlib):
#-------------------------------------------------------------------------------------------
 # define R user libray path
 if not os.path.exists(psicsRlib):
  os.makedirs(psicsRlib)

 if not os.environ.get("R_LIBS_USER"):
  os.environ.setdefault("R_LIBS_USER",psicsRlib)
 else:
  os.environ["R_LIBS_USER"] = psicsRlib

 # define psics package source code os path
 ref = resources.files('psics')/'BoxQP.mex'
 with resources.as_file(ref) as path:
  vpath = path
 path = utils.defpath(str(vpath))

 #-----------------------------------------------------------------------------------
 # setup covsel code
 #-----------------------------------------------------------------------------------
 print("psics: setup covsel code")
 if dwnCovsel:
  #-----------------------------------------------------------------------------------
  # download covsel code
  #-----------------------------------------------------------------------------------

  # https://www.di.ens.fr/~aspremon/ZIP/COVSEL.zip
  url = (urlCovsel)
  filename = "COVSEL.zip"
  urlretrieve(url, filename)

  # unzip covsel.zip file in directory covsel
  with zipfile.ZipFile("COVSEL.zip", 'r') as zip_ref:
   zip_ref.extractall("covsel")

  # workarounds to fix errors in original code and changes to compile BoXQP in octave

  # 1)  Change line #include <matrix.h> to #include <cblas.h> in file BoxQP.h
  modify("covsel/BoxQP/BoxQP.h","<matrix.h>","<cblas.h>")

  # 2) Fix code error in BoxQP_mex.c replace #include "BoxQp.h" with #include "BoxQP.h"
  modify("covsel/BoxQP/BoxQP_mex.c","BoxQp.h","BoxQP.h")

  # 3) Fix code error in utils.c, drop win32 shell value 
  modify("covsel/BoxQP/utils.c","#endif win32","#endif")

  # 4) compile BoXQP in octave
  pth = "covsel/BoxQP/"
  cp = subprocess.run(["mkoctfile","--mex",pth+"BoxQP.c",pth+"BoxQP_mex.c",
      pth+"utils.c"])
  if cp.returncode !=0:
   print("Error compiling: mkoctfile --mex BoxQP.c BoxQP_mex.c utils.c")
   quit()

  # 5)copy matlab covsel program to psics package code 
  cp = subprocess.run(["cp","covsel/spmlcdvec.m",path+"/."])
  if cp.returncode !=0:
   print("Error copying spmlcdvec.m to psics/code directory")
   quit()

  # 6) remove covsel downloaded and extracted files
  cp = subprocess.run(["rm","-r","covsel"])
  if cp.returncode !=0:
   print("Error deleting temp directory covsel")
   quit()
  cp = subprocess.run(["rm","COVSEL.zip"])
  if cp.returncode !=0:
   print("Error deleting COVSEL.zip file")
   quit()



 #-----------------------------------------------------------------------------------
 # use covsel code included in psics
 #-----------------------------------------------------------------------------------

 else:
  # 5') compile BoxQp using package source files previously downloaded from d'Aspremont site

  cp = subprocess.run(["mkoctfile","--mex",path+"/BoxQP.c",path+"/BoxQP_mex.c",path+"/utils.c"
                      ])
  if cp.returncode !=0:
   print("Error compiling previousily downloaded: mkoctfile --mex BoxQP.c BoxQP_mex.c utils.c")
   quit()

 #  move BoxQP.mex to package code
 cp = subprocess.run(["mv","BoxQP.mex",path+"/."])
 if cp.returncode !=0:
  print("Error moving BoxQP.mex psics/code directory")
  quit()



 #--------------------------------------------------------------------------------------
 # setup covpath code
 #--------------------------------------------------------------------------------------

 if dwnCovpath:
  # download covpath code
  print("psics: setup covpath code")
  url = (urlCovpath)
  filename = "Covpath_1.0.tar.gz"
  urlretrieve(url, filename)

  # untar Covpath_1.0.tar.gz
  tar = tarfile.open(filename, "r:gz")
  tar.extractall()
  tar.close()

  # 1) Change line to include useDynLib clause
  modify("Covpath/NAMESPACE","# Default NAMESPACE created by R","useDynLib(Covpath,bcdCorrector)")

  # 2) Install Covpath package in R
  cp = subprocess.run(["R","CMD","INSTALL","Covpath"])
  if cp.returncode !=0:
   print("Error installing: R CMD INSTALL Covpath")
   quit()

  # 3) remove Covpath downloaded and extracted files
  cp = subprocess.run(["rm","-r","Covpath"])
  if cp.returncode !=0:
   print("Error deleting temp directory Covpath")
   quit()
  cp = subprocess.run(["rm",filename])
  if cp.returncode !=0:
   print("Error deleting Covpath_1.0.tar.gz file")
   quit()

 # setup previously corrected Covpath code included in psics package 
 else:
  #  Install local Covpath package in R
  cp = subprocess.run(["R","CMD","INSTALL",path+"/Covpath"])
  if cp.returncode !=0:
   print("Error installing local: R CMD INSTALL Covpath")
   quit()

 #-----------------------------------------------------------------------------------
 # setup glasso R package downloading from cran project
 #-----------------------------------------------------------------------------------
 # always download glasso code from R cran project
 print("psics: setup glasso code")
 installstr = "install.packages('glasso',repos='"+urlGlasso+"')"
 cp = subprocess.run(["Rscript","-e",installstr])
 if cp.returncode != 0:
  print("Error installing glasso")
  quit()

 #-----------------------------------------------------------------------------------
 # setup alm code
 #-----------------------------------------------------------------------------------
 print("psics: setup alm code")
 if dwnAlm:
  #-----------------------------------------------------------------------------------
  # download  alm code
  #-----------------------------------------------------------------------------------
  url = (urlAlm)
  filename = "master.zip"
  urlretrieve(url, filename)

  # unzip master.zip file in directory alm
  with zipfile.ZipFile("master.zip", 'r') as zip_ref:
   zip_ref.extractall("alm")
  
  # copy matlab ALM program to psics package code 
  cp = subprocess.run(["cp","alm/SICS-master/SICS_ALM/SICS_ALM.m",path+"/."])
  if cp.returncode !=0:
   print("Error copying SICS_ALM.m to psics/code directory")
   quit()  

# remove alm downloaded and extracted files
  cp = subprocess.run(["rm","-r","alm"])
  if cp.returncode !=0:
   print("Error deleting temp directory alm")
   quit()
  cp = subprocess.run(["rm","master.zip"])
  if cp.returncode !=0:
   print("Error deleting master.zip file")
   quit()
