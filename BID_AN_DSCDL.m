clear,clc;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
addpath('GatingNetwork');
GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_MeanImage\';
GT_fpath = fullfile(GT_Original_image_dir, '*.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_NoisyImage\';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

%% load parameters and dictionary
load Data/params_flexible.mat par param;
load Data/Dict_DSCDL_backup_flexible_20160724T132738.mat Dict;
load Data/EMGM_8x8_100_knnNI2BS500Train_20160722T082406.mat;
par.cls_num = 100;
par.nInnerLoop = 3;

PSNR = [];
SSIM = [];
CCPSNR = [];
CCSSIM = [];
for i = 1 : im_num
    IMin = im2double(imread(fullfile(TT_Original_image_dir,TT_im_dir(i).name) ));
    IM_GT = im2double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
    S = regexp(TT_im_dir(i).name, '\.', 'split');
    IMname = S{1};
    [h,w,ch] = size(IMin);
    CCPSNR = [CCPSNR csnr( IMin*255,IM_GT*255, 0, 0 )];
    CCSSIM = [CCSSIM cal_ssim( IMin*255, IM_GT*255, 0, 0 )];
    fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(end), CCSSIM(end));
    % color or gray image
    if ch==1
        IMin_y = IMin;
    else
        % change color space, work on illuminance only
        IMin_ycbcr = rgb2ycbcr(IMin);
        IMin_y = IMin_ycbcr(:, :, 1);
        IMin_cb = IMin_ycbcr(:, :, 2);
        IMin_cr = IMin_ycbcr(:, :, 3);
    end
    %%
    nOuterLoop = 1;
    Continue = true;
    while Continue
        fprintf('Iter: %d \n', nOuterLoop);
        IMout_y = bscdl_BID_full(IMin_y,model,Dict,par,param);
        % Noise Level Estimation
        nSig = NoiseLevel(IMout_y*255);
        fprintf('The noise level is %2.4f.\n',nSig);
        if nSig < 0.0001 || nOuterLoop >= 5
            Continue = false;
        else
            nOuterLoop = nOuterLoop + 1;
            IMin_y = IMout_y;
        end
    end
    if ch==1
        IMout = IMout_y;
    else
        IMout_ycbcr = zeros(size(IMin));
        IMout_ycbcr(:, :, 1) = IMout_y;
        IMout_ycbcr(:, :, 2) = IMin_cb;
        IMout_ycbcr(:, :, 3) = IMin_cr;
        IMout = ycbcr2rgb(IMout_ycbcr);
    end
    %% output
    PSNR = [PSNR csnr( IMout*255, IM_GT*255, 0, 0 )];
    SSIM = [SSIM cal_ssim( IMout*255, IM_GT*255, 0, 0 )];
    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', PSNR(end), SSIM(end));
    %% output
    imwrite(IMout, ['../cc_Results/Real_DSCDL/DSCDL_BID_AN_NLM_' IMname '.png']);
end
mPSNR = mean(PSNR);
mSSIM = mean(SSIM);
mCCPSNR = mean(CCPSNR);
mCCSSIM = mean(CCSSIM);
save(['C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\DSCDL_BID_AN_NLM.mat'],'PSNR','mPSNR','SSIM','mSSIM','CCPSNR','mCCPSNR','CCSSIM','mCCSSIM');
