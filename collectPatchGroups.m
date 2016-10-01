clear;
addpath('Data');
addpath('Utilities');
% TrainingClean = '../TrainingData/ycbcrDenoised/';
TrainingNoisy = '../TrainingData/RGBNoisy/';
% TrainingClean = '../TrainingData/ycbcrDenoised/';
TrainingClean = '../TrainingData/RGBBSDS500train/';

patch_size = 8; 
Patch_Channel = 1;
num_patch_N = 200000;
num_patch_C = 5*num_patch_N;
R_thresh = 0.05;
cls_num = 64;

% Parameters Setting
par.nlsp = 10;
par.step               =    2;
par.Patch_Channel = Patch_Channel;
par.patch_size                =    patch_size;
par.S         =  2 * par.patch_size - 1;
par.rho = 5e-2;
par.lambda1         =       0.01;
par.lambda2         =       0.001;
par.mu              =       0.01;
par.sqrtmu          =       sqrt(par.mu);
par.nu              =       0.1;
par.nIter           =       100;
par.epsilon         =       5e-3;
par.t0              =       5;
par.K = 256;
par.L = par.patch_size * par.patch_size;
par.cls_num = cls_num;
param.K = par.K;
param.lambda = par.lambda1;
param.iter=300; 
param.L = par.patch_size * par.patch_size;
save Data/params_gray_PG.mat par param;

[XN, XC] = rnd_smp_PG_kNN(TrainingNoisy, TrainingClean, num_patch_N, num_patch_C,  par);

[model, llh, cls_idx]  =  empggm(XC, cls_num,par.nlsp); 
for c = 1 : cls_num
    idx = find(cls_idx == c);
%     if (length(idx) >  10000)
%         select_idx = randperm(length(idx));
%         idx = idx(select_idx(1:10000));
%     end
    Xn{c} = XN(:, idx);
    Xc{c} = XC(:, idx);
end
% % add smooth patches
% Xn{cls_num+1} = XN0; 
% Xc{cls_num+1} = XC0;
% [s_idx, seg]    =  Proc_cls_idx( cls_idx );
% cls_num = par.cls_num+1;
% model.means(:,cls_num) = mean(XC0,2);
% model.covs(:,:,cls_num) = cov(XC0');
% length0 = size(XC0,2)/par.nlsp;
% model.mixweights = [model.mixweights length0/(length0 + length(cls_idx))]/(sum(model.mixweights) + length0/(length0 + length(cls_idx)));
% model.nmodels = model.nmodels + 1;
% clear XN XC XN0 XC0;
% save model
GMM_model = ['Data/GMM_PG_' num2str(par.nlsp) '_' num2str(patch_size) 'x' num2str(patch_size) '_' num2str(cls_num) '_' datestr(now, 30) '.mat'];
save(GMM_model, 'model', 'Xn','Xc','cls_num');  