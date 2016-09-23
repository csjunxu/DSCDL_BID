%       3) when param.mode=2
%         min_{alpha} 0.5||x-Dalpha||_2^2 + lambda||alpha||_1 +0.5 lambda2||alpha||_2^2
% key parameter for performance:
% par.sqrtmu     param.lambda       param.lambda2

clear;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_MeanImage\';
GT_fpath = fullfile(GT_Original_image_dir, '*.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_NoisyImage\';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

%% load parameters and dictionary
load Data/params.mat par param;
load Data/Dict_DSCDL_backup_flexible_20160724T132738.mat Dict;
load Data/EMGM_8x8_100_knnNI2BS500Train_20160722T082406.mat;
par.cls_num = 100;
for lambda = [0.01 0.05 0.1:0.1:1]
    param.lambda = lambda;
    for lambda2 = [0.0001 0.001 0.01 0.1 0.5 1]
        param.lambda2 = lambda2;
        for sqrtmu = 0.1
            par.sqrtmu = sqrtmu;
            for nInnerLoop = [3 4]
                par.nInnerLoop = nInnerLoop;
                load ImageIndex.mat;
                PSNR = [];
                SSIM = [];
                CCPSNR = [];
                CCSSIM = [];
                for i = 1:im_num
                    IMmean = im2double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
                    IMin = im2double(imread(fullfile(TT_Original_image_dir, TT_im_dir(i).name)));
                    S = regexp(TT_im_dir(i).name, '\.', 'split');
                    IMname = S{1};
                    fprintf('%s : \n',IMname);
                    CCPSNR = [CCPSNR csnr( IMin*255,IMmean*255, 0, 0 )];
                    CCSSIM = [CCSSIM cal_ssim( IMin*255, IMmean*255, 0, 0 )];
                    fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n',CCPSNR(end), CCSSIM(end));
                    [h,w,ch] = size(IMin);
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
                        % fprintf('Iter: %d \n', nOuterLoop);
                        IMout_y = bscdl_BID_full(IMin_y,model,Dict,par,param);
                        % Noise Level Estimation
                        nSig = NoiseLevel(IMout_y*255);
                        fprintf('The noise level is %2.4f.\n',nSig);
                        if nSig < 0.0001 && nOuterLoop >= 1
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
                    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', csnr( IMout*255,IMmean*255, 0, 0 ), cal_ssim( IMout*255, IMmean*255, 0, 0 ));
                    %% output
                    imwrite(IMout, ['../cc_Results/Real_ours/DSCDL_BID_AN_CCNoise_Conv_' IMname '_' num2str(sqrtmu) '_lambda_' num2str(lambda) '_lambda2_' num2str(lambda2) '.png']);
                    PSNR =[PSNR csnr( IMout*255,IMmean*255, 0, 0 )];
                    SSIM = [SSIM cal_ssim( IMout*255, IMmean*255, 0, 0 )];
                    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', PSNR(end), SSIM(end));
                end
                mPSNR = mean(PSNR);
                mSSIM = mean(SSIM);
                savename = ['../cc_Results/DSCDL_BID_AN_CCNoise_nInLoop' num2str(nInnerLoop) '_sqrtmu_' num2str(sqrtmu) '_lambda_' num2str(lambda) '_lambda2_' num2str(lambda2) '.mat'];
                save(savename, 'mPSNR', 'mSSIM', 'PSNR', 'SSIM');
            end
        end
    end
end