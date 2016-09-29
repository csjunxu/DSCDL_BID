clear;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
addpath('GatingNetwork');
Original_image_dir = './TestingImages/';
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

%% load parameters and dictionary
load Data/params.mat par param;
load Data/DSCDL_BID_Dict_6x6x3_ADPU_backup_nup0_20160927T221836.mat Dict;
load Data/EMGM_6x6x3_100_20160927T093144.mat;
par.cls_num = 100;
par.nInnerLoop = 3;

for i = 1 : im_num
    IMin=im2double(imread(fullfile(Original_image_dir, im_dir(i).name)));
    S = regexp(im_dir(i).name, '\.', 'split');
    IMname = S{1};
    [h,w,ch] = size(IMin);
    %%
    nOuterLoop = 1;
    Continue = true;
    while Continue
        fprintf('Iter: %d \n', nOuterLoop);
        IMout = bscdl_BID_RGB(IMin,model,Dict,par,param);
        % Noise Level Estimation
        nSig = NoiseEstimation(im2uint8(IMout),par.win);
        fprintf('The noise level is %2.4f.\n',nSig);
        if nSig < 0.0001 || nOuterLoop >= 10
            Continue = false;
        else
            nOuterLoop = nOuterLoop + 1;
            IMin = IMout;
        end
    end
    % output
    imwrite(IMout, ['./DSCDL_BID_AN_ADPU/ADPU_nup0_DSCDL_BID_AN_' IMname '.png']);
    fprintf('The Loops of these parts are  %s. \n', num2str(nOuterLoop));
end