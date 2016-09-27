clear;clc;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
% addpath('SPAMS/release/mkl64');

load Data/params.mat;
load Data/EMGM_6x6x3_100_20160927T093144.mat;
% Parameters Setting
par.rho = 0.05;
par.lambda1         =       0.01;
par.lambda2         =       0.01;
par.mu              =       0.01;
par.sqrtmu          =       sqrt(par.mu);
par.nu              =       0.1;
par.nup              =      0;
par.win                =    6;
par.channel = 3;
% fixed parameters
par.epsilon         =        5e-3; 
par.cls_num            =    cls_num;
par.step               =    2;
par.nIter           =       100;
par.t0              =       5;
par.K               =       256;
par.L               =       par.win * par.win * par.channel;
% parameters for mexLasso
param.K = par.K;
param.iter=300;
param.lambda = par.lambda1;
param.lambda2 = par.lambda2;
param.L = par.L;
flag_initial_done = 0;
paramsname = sprintf('Data/params.mat');
save(paramsname,'par','param');

load Data/EMGM_6x6x3_100_20160927T093144.mat;
% Initiate Dictionary
Dini = [];
for i = 1 : par.cls_num
    XN_t = double(Xn{i});
    XC_t = double(Xc{i});
    XN_t = XN_t - repmat(mean(XN_t), [size(XN_t,1) 1]);
    XC_t = XC_t - repmat(mean(XC_t), [size(XC_t,1) 1]);
    fprintf('Double Semi-Coupled dictionary learning: Cluster: %d\n', i);
    D = mexTrainDL([XN_t;XC_t], param);
    Dini{i} = D;
    Dict_BID_Initial = sprintf('Data/Dict_DSCDL_BID_Dict_6x6x3_ADPU_Initial_nup0_%s.mat', datestr(now, 30));
    save(Dict_BID_Initial,'Dini');
    D = Dini{i};
    clear Dini;
    Dn = D(1:size(XN_t,1),:);
    Dc = D(size(XN_t,1)+1:end,:);
    Wn = eye(size(Dn, 2));
    Wc = eye(size(Dc, 2));
    Alphac = mexLasso([XN_t;XC_t], D, param);
    Alphan = Alphac;
    clear D;
    fprintf('Double Semi-Coupled dictionary learning: Cluster: %d\n', i);
    [Alphac, Alphan, XC_t, XN_t, Dc, Dn, Wc, Wn, Uc, Un, Pn, f] = ADPU_Double_Semi_Coupled_DL(Alphac, Alphan, XC_t, XN_t, Dc, Dn, Wc, Wn, par);
    Dict.DC{i} = Dc;
    Dict.DN{i} = Dn;
    Dict.WC{i} = Wc;
    Dict.WN{i} = Wn;
    Dict.UC{i} = Uc;
    Dict.UN{i} = Un;
    Dict.PN{i} = Pn;
    Dict.f{i} = f;
    Dict_BID_backup = sprintf('Data/DSCDL_BID_Dict_6x6x3_ADPU_backup_nup0_%s.mat',datestr(now, 30));
    save(Dict_BID_backup,'Dict');
end