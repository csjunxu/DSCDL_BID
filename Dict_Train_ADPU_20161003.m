clear;clc;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
% addpath('SPAMS/release/mkl64');

load Data/params_gray_PG.mat;
load Data/GMM_PG_3_10_8x8_64_20161003T094301.mat;
% Parameters Setting
par.rho = 0.05;
par.lambda1         =       0.01;
par.lambda2         =       0.001;
par.mu              =       0.01;
par.sqrtmu          =       sqrt(par.mu);
par.nu              =       0.1;
par.nup              =      0;
par.epsilon         =        5e-3;
par.cls_num            =    cls_num;
par.step               =    2;
par.win                =    par.patch_size;
par.nIter           =       100;
par.t0              =       5;
par.K               =       256;
par.L               =       par.win * par.win;
param.K = par.K;
param.iter=300;
param.lambda = par.lambda1;
param.lambda2 = par.lambda2;
param.L = par.win * par.win;
flag_initial_done = 0;
paramsname = sprintf('Data/params_gray_PG.mat');
save(paramsname,'par','param');

load Data/GMM_PG_10_8x8_64_20160930T171410.mat;
% Initiate Dictionary
Dini = cell(1,size(model,1));
for j = 1:size(model,1)
    for i = 1 : par.cls_num
        XN_t = double(Xn{j,i});
        XC_t = double(Xc{j,i});
        XN_t = XN_t - repmat(mean(XN_t), [par.win^2 1]);
        XC_t = XC_t - repmat(mean(XC_t), [par.win^2 1]);
        fprintf('Double Semi-Coupled dictionary learning: Cluster: %d\n', i);
        D = mexTrainDL([XN_t;XC_t], param);
        Dini{j,i} = D;
        Dict_BID_Initial = sprintf('Data/DSCDL_Dict_PG_3_10_8x8_64_BID_Initial_%s.mat', datestr(now, 30));
        save(Dict_BID_Initial,'Dini');
        D = Dini{j,i};
        clear Dini;
        Dn = D(1:par.win * par.win,:);
        Dc = D(par.win * par.win+1:end,:);
        Wn = eye(size(Dn, 2));
        Wc = eye(size(Dc, 2));
        Alphac = mexLasso([XN_t;XC_t], D, param);
        Alphan = Alphac;
        clear D;
        fprintf('Double Semi-Coupled dictionary learning: Cluster: %d\n', i);
        [Alphac, Alphan, XC_t, XN_t, Dc, Dn, Wc, Wn, Uc, Un, Pn, f] = ADPU_Double_Semi_Coupled_DL(Alphac, Alphan, XC_t, XN_t, Dc, Dn, Wc, Wn, par);
        Dict.DC{j,i} = Dc;
        Dict.DN{j,i} = Dn;
        Dict.WC{j,i} = Wc;
        Dict.WN{j,i} = Wn;
        Dict.UC{j,i} = Uc;
        Dict.UN{j,i} = Un;
        Dict.PN{j,i} = Pn;
        Dict.f{j,i} = f;
        Dict_BID_backup = sprintf('Data/DSCDL_ADPU_Dict_PG_3_10_8x8_64_BID_backup_%s.mat',datestr(now, 30));
        save(Dict_BID_backup,'Dict');
    end
end