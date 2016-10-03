clear;
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
par.nInnerLoop = 1;

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
    fprintf('%s: \n',TT_im_dir(i).name);
    CCPSNR = [CCPSNR csnr( IMin*255,IM_GT*255, 0, 0 )];
    CCSSIM = [CCSSIM cal_ssim( IMin*255, IM_GT*255, 0, 0 )];
    fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(end), CCSSIM(end));
    %% 3 channels
    IMout = zeros(size(IMin));
    for cc = 1:ch
        IMin_cc = IMin(:,:,cc);
        IM_GT_cc = IM_GT(:,:,cc);
        fprintf('Channel %d: The initial PSNR = %2.4f, SSIM = %2.4f. \n', cc, csnr( IMin_cc*255,IM_GT_cc*255, 0, 0 ), cal_ssim( IMin_cc*255, IM_GT_cc*255, 0, 0 ));
        IMout_cc = bscdl_BID_full(IMin_cc,model,Dict,par,param);
        fprintf('Channel %d: The medium PSNR = %2.4f, SSIM = %2.4f. \n', cc, csnr( IMout_cc*255, IM_GT_cc*255, 0, 0 ), cal_ssim( IMout_cc*255, IM_GT_cc*255, 0, 0 ));
        nSig = 10;
        [PGPDpar, model]  =  Parameters_Setting( nSig );
        PGPDpar.I = IM_GT_cc;
        PGPDpar.nim = IMout_cc;
        PGPDpar.IteNum = 1;
        [IMout_cc,PGPDpar]  =  PGPD_Denoising(PGPDpar,model);
        fprintf('Channel %d: The final PSNR = %2.4f, SSIM = %2.4f. \n', cc, csnr( IMout_cc*255, IM_GT_cc*255, 0, 0 ), cal_ssim( IMout_cc*255, IM_GT_cc*255, 0, 0 ));
        IMout(:,:,cc) = IMout_cc;
    end
    %% output
    PSNR = [PSNR csnr( IMout*255, IM_GT*255, 0, 0 )];
    SSIM = [SSIM cal_ssim( IMout*255, IM_GT*255, 0, 0 )];
    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', csnr( IMout*255, IM_GT*255, 0, 0 ), cal_ssim( IMout*255, IM_GT*255, 0, 0 ));
    %% output
    imwrite(IMout, ['../cc_Results/Real_DSCDL/DSCDL_P_PGPD_3Chs_BID_InLoop_' num2str(par.nInnerLoop) '_' IMname '.png']);
end
mPSNR = mean(PSNR);
mSSIM = mean(SSIM);
mCCPSNR = mean(CCPSNR);
mCCSSIM = mean(CCSSIM);
save(['C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\DSCDL_P_PGPD_3Chs_BID_InLoop_1.mat'],'PSNR','mPSNR','SSIM','mSSIM','CCPSNR','mCCPSNR','CCSSIM','mCCSSIM');
