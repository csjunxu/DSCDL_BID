function [XN, XC, XN0, XC0] = rnd_smp_PG_couple(TrainingNoisy, TrainingClean, nSig, par, num_patch, R_thresh)

Nim_path = fullfile(TrainingNoisy,'*.tif');
Cim_path = fullfile(TrainingClean,'*.tif');

Nim_dir = dir(Nim_path);
Cim_dir = dir(Cim_path);

Nim_num = length(Nim_dir);
Cim_num = length(Cim_dir);
nper_img = zeros(1, Cim_num);
for ii = 1:Cim_num
    Cim = im2double(imread(fullfile(TrainingClean, Cim_dir(ii).name)));
    [h,w] = size(Cim);
    if h >= 1000
        randh = randi(h-1000);
        Cim = Cim(randh+1:randh+1000,:,:);
    end
    if w >= 1000
        randw = randi(w-1000);
        Cim = Cim(:,randw+1:randw+1000,:);
    end
    nper_img(ii) = numel(Cim);
end

nper_img = floor(nper_img*num_patch/sum(nper_img));

XN = [];
XC = [];
XN0 = [];
XC0 = [];
for ii = 1:Cim_num
    warning off;
    patch_num = nper_img(ii);
    Cim = im2double(imread(fullfile(TrainingClean, Cim_dir(ii).name)));
    % generate noisy image
    randn('seed',0);
    Nim =   Cim + nSig/255*randn(size(Cim));
    fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', csnr( Nim*255, Cim*255, 0, 0 ),cal_ssim( Nim*255, Cim*255, 0, 0 ));

    [h,w] = size(Nim);
    if h >= 1000
        randh = randi(h-1000);
        Nim = Nim(randh+1:randh+1000,:,:);
        Cim = Cim(randh+1:randh+1000,:,:);
    end
    if w >= 1000
        randw = randi(w-1000);
        Nim = Nim(:,randw+1:randw+1000,:);
        Cim = Cim(:,randw+1:randw+1000,:);
    end
    
    [NPG, CPG,NPG0, CPG0] = sample_PG_couple(Nim, Cim, par, patch_num, R_thresh);
    XN = [XN, NPG];
    XC = [XC, CPG]; 
    XN0 = [XN0, NPG0];
    XC0 = [XC0, CPG0]; 
end
num_patch = size(XN,2);

patch_path = ['Data/rnd_PG_' num2str(par.patch_size) 'x' num2str(par.patch_size) '_' num2str(num_patch) '_' num2str(R_thresh)  '_' datestr(now, 30) '.mat'];
save(patch_path, 'XN', 'XC','XN0','XC0');