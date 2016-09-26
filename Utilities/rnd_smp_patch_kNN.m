function [XN, XC, patch_size, ch] = rnd_smp_patch_kNN(TrainingNoisy, TrainingClean, num_patch_N, num_patch_C, Patch_Channel, R_thresh)
warning off;
Nim_path = fullfile(TrainingNoisy,'*.jpg');
Cim_path = fullfile(TrainingClean,'*.jpg');

Nim_dir = dir(Nim_path);
Cim_dir = dir(Cim_path);

Nim_num = length(Nim_dir);
Cim_num = length(Cim_dir);

% noisy patches per image
nper_img_N = zeros(1, Nim_num);
for ii = 1:Nim_num
    Nim = im2double(imread(fullfile(TrainingNoisy, Nim_dir(ii).name)));
    [h,w,ch] = size(Nim);
    if h >= 1000
        randh = randi(h-1000);
        Nim = Nim(randh+1:randh+1000,:,:);
    end
    if w >= 1000
        randw = randi(w-1000);
        Nim = Nim(:,randw+1:randw+1000,:);
    end
    nper_img_N(ii) = numel(Nim);
end
nper_img_N = floor(nper_img_N*num_patch_N/sum(nper_img_N));

% clean patches per image
nper_img_C = zeros(1, Cim_num);
for ii = 1:Cim_num
    Cim = im2double(imread(fullfile(TrainingClean, Cim_dir(ii).name)));
    [h,w,ch] = size(Cim);
    if h >= 1000
        randh = randi(h-1000);
        Cim = Cim(randh+1:randh+1000,:,:);
    end
    if w >= 1000
        randw = randi(w-1000);
        Cim = Cim(:,randw+1:randw+1000,:);
    end
    nper_img_C(ii) = numel(Cim);
end
nper_img_C = floor(nper_img_C*num_patch_C/sum(nper_img_C));

% extract clean patches
XCa = [];
for ii = 1:Cim_num
    patch_num = nper_img_C(ii);
    Cim = im2double(imread(fullfile(TrainingClean, Cim_dir(ii).name)));
    [h,w,ch] = size(Cim);
    if h >= 1000
        randh = randi(h-1000);
        Cim = Cim(randh+1:randh+1000,:,:);
    end
    if w >= 1000
        randw = randi(w-1000);
        Cim = Cim(:,randw+1:randw+1000,:);
    end
    [C, patch_size, ch] = sample_patches(Cim, patch_num, Patch_Channel, R_thresh);
    XCa = [XCa, C];
end

% extract noisy patches
XN = [];
XC = [];
for ii = 1:Nim_num
    patch_num = nper_img_N(ii);
    Nim = im2double(imread(fullfile(TrainingNoisy, Nim_dir(ii).name)));
    [h,w] = size(Nim);
    if h >= 1000
        randh = randi(h-1000);
        Nim = Nim(randh+1:randh+1000,:,:);
    end
    if w >= 1000
        randw = randi(w-1000);
        Nim = Nim(:,randw+1:randw+1000,:);
    end
    [N, patch_size, ch] = sample_patches(Nim, patch_num, Patch_Channel, R_thresh);
    XN = [XN, N];
    % given noisy patches, search corresponding clean ones via k-NN
    IDX = knnsearch(N', XCa');
    XC = [XC, XCa(:, IDX)];
    
%     % given noisy patches, search corresponding clean ones via pdist2
%     D = pdist2(N', XCa');
%     [~, IDX2] = sort(D,2);
%     kneighbours = IDX2(:,1);
%     XC2 = [XC2, XCa(:, kneighbours)];
end
% final results
num_patch = size(XN,2);
patch_path = ['Data/rnd_patches_' num2str(patch_size) 'x' num2str(patch_size) 'x' num2str(ch) '_' num2str(num_patch) '_' num2str(R_thresh) '_' datestr(now, 30) '.mat'];
save(patch_path, 'XN', 'XC');  