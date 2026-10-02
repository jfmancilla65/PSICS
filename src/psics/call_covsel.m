% psics python package
% file: call_covsel.m
% calls covsel algorithm
%
% read parameters  from command line
arg_list   = argv ();
sigma_file = arg_list{1};
rho        = str2double(arg_list{2});
maxiter    = str2num(arg_list{3});
prec       = str2double(arg_list{4});
maxnest    = str2num(arg_list{5});
algo       = arg_list{6};
packdir    = arg_list{7};

% add python package dir to path
addpath(packdir)

% read covariance matrix file
sigma = csvread(sigma_file);

% call COVSEL code
[X,U,gvals,cputimes] = spmlcdvec(sigma,rho,maxiter,prec,maxnest,algo);

% output file results
csvwrite('_covselXsol.dat',X);
csvwrite('_covselUsol.dat',U);
csvwrite('_covselgvals.dat',gvals);
csvwrite('_covselcputimes.dat',cputimes);


