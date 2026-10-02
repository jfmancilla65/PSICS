# PSICS

## A Python toolbox for numerically solving the sparse inverse covariance selection. 

## Introduction

This Python module includes integrated access to native code for **COVSEL**, **COVPATH**, **GLASSO**, and **ALM** algorithms that were previously published. 

## PSICS description

PSICS is built using the Python programming language.  Consists of several functions that integrate the call to each of the sparse inverse covariance selection (SICS) algorithms, using the standard library **subprocess**.

This approach makes it possible to call different software scripts or programs written in **Matlab** and **R**, the two main platforms used to build original algorithms.

PSCIS installs the native code algorithms **COVSEL**, **COVPATH**, **GLASSO** and **ALM**. Originally, **COVSEL** and ALM runned in a **Matlab** environment. In this package, **Matlab** was replaced with the Octave platform. 

## PSICS general requeriments

Next are listed the software general requirements to install and use PSICS. It is important to note that this Python package was built and tested on a Debian 12 linux system.

- Rscript (R) >= 4.2.2
- GNU Octave, >=  7.3.0 
- gcc and GNU Fortran >= 12.2.0
- blas and cblas libraries (included in Octave and R installation)
- Python >= 3.10.
- pip >= 23.0.1.

An active connection to the internet is required in order to download software.

## Directions to install  and configure PSICS

### Install general requirements 

It is necessary to install the R, Octave, gcc and gFORTRAN compilers. This can be accomplished using the following commands in the terminal of a Debian-based Linux system.  

#### The gcc and gfortran compiler
```
$ sudo apt-get install gcc gfortran
```

#### R software
```
$ sudo apt-get install r-base r-base-core
```

#### Octave software
```
$ sudo apt-get install octave octave-dev

```

### Create a Python virtual environment

In order to avoid any conflict with the user software configuration, a virtual environment must be created and activated to install and test PSICS.

```
$ python3 -m venv ~/.venv
$ source ~/.venv/bin/activate
```


### PSICS installation

PSICS can be installed from pypi site.

```
$ pip install psics
```


##### Configure the PSICS package

In this step  the algorithms COVPATH, COVSEL, ALM and GLASSO native programs are configured.

``` 
$ python3 
>>> import psics.config as config
>>> config.setup()
>>> ...
>>> quit()
```



## Testing PSICS package


To verify the usability of the PSICS code, a testing  script is included.

In this script, the following tasks are executed:

* Create a synthetic covariance matrix with a previously known structure

* Keep the original synthetic inverse as a reference to compare the results obtained from SICS algorithms

* Call ALM, COVSEL, COVPATH, GLASSO and GLASSOPATH algorithms

* Calculate some key indicators for the results obtained in each case

* Show results resume in a text table and some graphical outputs where the real inverse covariance is compared to the inverse covariance calculated

```
$ python3 
>>> import psics.test as test
>>> test.run()
>>> ...
>>> quit()
```

Further Python scripts to use PSICS can be coded using this program test as reference.


## Workarounds applied 

Because of some technical issues, it was required to apply some workarounds to COVPATH and COVSEL code.


1. COVPATH is not available to be installed as an R package directly from the CRAN Project site. It seems to be a version compatibility problem in the use of the useDynLib clause in NAMESPACE configuration file. Internal PSICS setup.py function corrects this problem.


2. COVSEL original code included the Quadratic program solver BoxQP compiled as a Matlab mex file. Because PSICS uses Octave instead of Matlab, setup.py function replaces the Matlab mex by compiling BoxQP in the Octave mkoctfile standard.


3. Some other minor errors in the COVSEL code source are also corrected automatically by the setup.py function.


## User's guide

Detailed guide of use can be found in file USERSGUIDE.md included in this Python package.

