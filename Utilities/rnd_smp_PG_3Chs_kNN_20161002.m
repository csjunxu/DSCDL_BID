function [XN, XC] = rnd_smp_PG_3Chs_kNN_20161002(TrainingNoisy, TrainingClean, num_patch_N, num_patch_C, par)

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
CPGall = {};
CPGallmean  = {};
% extract noisy PGs and corresponding clean PGs
XN = {};
XNmean = {};
XC = {};
XCmean = {};
for cc = 1:ch
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
        Cim = Cim(:,:,cc);
        [CPG, CPGmean] = sample_PatchGroups(Cim, patch_num, par);
        CPGall{cc} = [CPGall{cc}, CPG];
        CPGallmean{cc} = [CPGallmean{cc}, CPGmean];
    end
    for ii = 1:Nim_num
        patch_num = nper_img_N(ii);
        Nim = im2double(imread(fullfile(TrainingNoisy, Nim_dir(ii).name)));
        [h,w,ch] = size(Nim);
        if h >= 1000
            randh = randi(h-1000);
            Nim = Nim(randh+1:randh+1000,:,cc);
        end
        if w >= 1000
            randw = randi(w-1000);
            Nim = Nim(:,randw+1:randw+1000,cc);
        end
        [NPG, NPGmean] = sample_PatchGroups(Nim, patch_num, par);
        XN{cc} = [XN{cc}, NPG];
        XNmean{cc} = [XNmean{cc}, NPGmean];
        % given noisy patches, search corresponding clean ones via k-NN
        NP = NPG(:,1:par.nlsp:end);
        CP  = CPGall(:,1:par.nlsp:end);
        PIDX = knnsearch(NP', CP');
        PIDX = (PIDX-1)*par.nlsp+1;
        PGIDX = PIDX';
        for jj = 1:par.nlsp-1
            PGIDX = [PGIDX;PIDX'+jj];
        end
        PGIDX = PGIDX(:);
        XC{cc} = [XC{cc}, CPGall(:, PGIDX)];
        XCmean{cc} = [XCmean{cc}, CPGallmean(:, PGIDX)];
    end
    XN{cc} = XN{cc} - XNmean{cc};
    XC{cc} = XC{cc} - XCmean{cc};
    num_patch = size(XN,2);
end
patch_path = ['Data/rnd_PG_' num2str(par.patch_size) 'x' num2str(par.patch_size) '_' num2str(num_patch)  '_' datestr(now, 30) '.mat'];
save(patch_path, 'XN', 'XC');