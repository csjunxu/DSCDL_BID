clear;
addpath('Data');
addpath('Utilities');
TrainingNoisy = '../TrainingData/RGBNoisy/';
% TrainingClean = '../TrainingData/ycbcrDenoised/';
TrainingClean = '../TrainingData/RGBBSDS500train/';

Patch_Channel = 1;
num_patch_N = 100000;
num_patch_C = 5*num_patch_N;
R_thresh = 0.05;
cls_num = 100;

[XN, XC, patch_size, ch] = rnd_smp_patch_kNN(TrainingNoisy, TrainingClean, num_patch_N, num_patch_C, Patch_Channel, R_thresh);

run('Param_Setting.m');
% load TrainingImages/rnd_patches_8x8_866242_0.05_20160715T232001.mat

[model, cls_idx]  =  emgm(XC, cls_num); 
for c = 1 : cls_num
    idx = find(cls_idx == c);
    if (length(idx) >  100000)
         select_idx = randperm(length(idx));
        idx = idx(select_idx(1:100000));
    end
    Xn{c} = XN(:, idx);
    Xc{c} = XC(:, idx);
end
clear XN XC;
GMM_model = ['Data/EMGM_' num2str(patch_size) 'x' num2str(patch_size) 'x'  num2str(ch) '_' num2str(cls_num) '_' datestr(now, 30) '.mat'];
save(GMM_model, 'model', 'Xn','Xc','cls_num');  