% Parameters Setting
par.cls_num            =    100;
par.step               =    2;
par.win                =    8;
par.rho = 5e-2;
par.lambda1         =       0.01;
par.lambda2         =       0.001;
par.mu              =       0.01;
par.sqrtmu          =       sqrt(par.mu);
par.nu              =       0.1;
par.nIter           =       100;
par.epsilon         =       5e-3;
par.t0              =       5;
par.K               =       256;
par.L               =       par.win * par.win;
param.K = par.K;
param.lambda = par.lambda1;
param.iter=300; 
param.L = par.win * par.win;
flag_initial_done = 0;

save params.mat par param;