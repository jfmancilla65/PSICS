% psics python package
% file: call_alm.m
% calls alm algorithm passing parameters read from command line
%
% read parameters from octave command line 
arg_list      = argv ();
cov_file      = arg_list{1};
rho           = str2double(arg_list{2});
opts.mxitr    = str2num(arg_list{3}); 
opts.mu0      = str2double(arg_list{4}); 
opts.muf      = str2double(arg_list{5}); 
opts.rmu      = str2double(arg_list{6});
opts.tol_gap  = str2double(arg_list{7});  
opts.tol_frel = str2double(arg_list{8});   
opts.tol_Xrel = str2double(arg_list{9});  
opts.tol_Yrel = str2double(arg_list{10});  
opts.tol_pinf = str2double(arg_list{11});  
opts.numDG    = str2num(arg_list{12});
opts.record   = str2num(arg_list{13});
opts.sigma    = str2double(arg_list{14});
packdir = arg_list{15};

% add packdir to path
addpath(packdir)

% read covariance matrix file
B = csvread(cov_file);

% call ALM algorithm
out = SICS_ALM(B,rho,opts);

% get results
res = zeros(6,1);
res(1,1) = out.iter ; 
res(2,1) = out.pinf ;
res(3,1) = out.obj ;
res(4,1) = out.gapX ; 
res(5,1) = out.gapY ;
res(6,1) = out.gap;

% output file results
csvwrite('_almXsol.dat',out.X);
csvwrite('_almYsol.dat',out.Y);
csvwrite('_almres.dat',res);
