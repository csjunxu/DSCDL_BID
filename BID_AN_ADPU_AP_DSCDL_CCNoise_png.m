clear;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
Original_image_dir = 'C:/Users/csjunxu/Desktop/Projects/crosschannel_CVPR2016/ccnoise_denoising/';
fpath = fullfile(Original_image_dir, '*_d800_iso6400_3_real.png');
GT_fpath = fullfile(Original_image_dir, '*_mean.png');
CC_fpath = fullfile(Original_image_dir, '*_ours.png');
im_dir  = dir(fpath);
GT_im_dir  = dir(GT_fpath);
CC_im_dir  = dir(CC_fpath);
im_num = length(im_dir);

%% load parameters and dictionary
load Data/params.mat par param;
load Data/DSCDL_BID_Dict_ADPU_backup_nup0.5_20160819T042830.mat Dict;
load Data/EMGM_8x8_100_knnNI2BS500Train_20160722T082406.mat;
par.cls_num = 100;
par.nInnerLoop = 1;
%% the whole image or part
type = 'all';
% 'random';
% 'middle';
% 'all';
PSNR = zeros(1,im_num);
SSIM = zeros(1,im_num);
CCPSNR = zeros(1,im_num);
CCSSIM = zeros(1,im_num);
for nup = 0.2:0.1:2
    par.nup = nup;
    for i = 1 : im_num
        IMin = im2double(imread(fullfile(Original_image_dir, im_dir(i).name)));
        S = regexp(im_dir(i).name, '\.', 'split');
        IMname = S{1};
        IMmean = im2double(imread(fullfile(Original_image_dir, GT_im_dir(i).name)));
        IM_CC = im2double(imread(fullfile(Original_image_dir, CC_im_dir(i).name)));
        fprintf('%s : \n',IMname);
        fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n',csnr( IMin*255, IMmean*255, 0, 0 ),cal_ssim( IMin*255, IMmean*255, 0, 0 ));
        [h,w,ch] = size(IMin);
        hh = [0:1000:h,h];
        ww = [0:1000:w,w];
        fprintf('There are/is %d parts in the image %s.%s!\n',(length(hh)-1)*(length(ww)-1),IMname,S{2});
        num_part = 0;
        if strcmp(type, 'all')
            listh = 1 : length(hh)-1;
            listw = 1 : length(ww)-1;
            IMout = zeros(h,w,ch);
        elseif strcmp(type, 'random')
            listh = randi([1, length(hh)-1],1,1);
            listw = randi([1, length(ww)-1],1,1);
            IMout = zeros(1000,1000,3);
        elseif strcmp(type, 'middle')
            listh = floor(median(1:length(hh)-1));
            listw = floor(median(1:length(ww)-1));
            IMout = zeros(1000,1000,3);
        end
        %%
        for nh = listh
            for nw = listw
                num_part = num_part + 1;
                IMin_part = IMin(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:);
                IMmean_part = IMmean(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:);
                % color or gray image
                if ch==1
                    IMin_part_y = IMin_part;
                else
                    % change color space, work on illuminance only
                    IMin_part_ycbcr = rgb2ycbcr(IMin_part);
                    IMin_part_y = IMin_part_ycbcr(:, :, 1);
                    IMin_part_cb = IMin_part_ycbcr(:, :, 2);
                    IMin_part_cr = IMin_part_ycbcr(:, :, 3);
                end
                %%
                nOuterLoop = 1;
                Continue = true;
                while Continue
                    fprintf('Iter: %d \n', nOuterLoop);
                    %                 IMout_part_y = bscdl_BID_full(IMin_part_y,model,Dict,par,param);
                    [IMout_part_y, PN] = AP_bscdl_BID_full(IMin_part_y,model,Dict,par,param);
                    % Noise Level Estimation
                    nSig = NoiseLevel(IMout_part_y*255);
                    fprintf('The noise level is %2.4f.\n',nSig);
                    if nSig < 0.0001 || nOuterLoop >= 1
                        Continue = false;
                    else
                        nOuterLoop = nOuterLoop + 1;
                        IMin_part_y = IMout_part_y;
                    end
                    if ch==1
                        IMout_part = IMout_part_y;
                    else
                        IMout_part_ycbcr = zeros(size(IMin_part));
                        IMout_part_ycbcr(:, :, 1) = IMout_part_y;
                        IMout_part_ycbcr(:, :, 2) = IMin_part_cb;
                        IMout_part_ycbcr(:, :, 3) = IMin_part_cr;
                        IMout_part = ycbcr2rgb(IMout_part_ycbcr);
                    end
                    fprintf('This is the %d/%d part of the image %s.%s!\n',num_part,(length(hh)-1)*(length(ww)-1),IMname,S{2});
                    if strcmp(type, 'all')
                        IMout(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:) = IMout_part;
                        IMmean(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:) = IMmean_part;
                    elseif strcmp(type, 'random') || strcmp(type, 'middle')
                        IMout = IMout_part;
                        IMmean = IMmean_part;
                    end
                    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', csnr( IMout_part*255,IMmean_part*255, 0, 0 ), cal_ssim( IMout_part*255, IMmean_part*255, 0, 0 ));
                end
            end
        end
        %% output
        imwrite(IMout, ['./DSCDL_BID_AN_ADPU_AP_CCNoise/DSCDL_BID_AN_ADPU_AP_CCNoise_' type '_' IMname '_' num2str(nup) '.png']);
        PSNR(i) = csnr( IMout*255,IMmean*255, 0, 0 );
        SSIM(i) = cal_ssim( IMout*255, IMmean*255, 0, 0 );
        CCPSNR(i) = csnr( IM_CC*255,IMmean*255, 0, 0 );
        CCSSIM(i) = cal_ssim( IM_CC*255, IMmean*255, 0, 0 );
        fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', PSNR(i), SSIM(i));
        fprintf('The CC:  PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(i), CCSSIM(i));
    end
%     mPSNR = mean(PSNR);
%     mSSIM = mean(SSIM);
%     CCmPSNR = mean(CCPSNR);
%     CCmSSIM = mean(CCSSIM);
%     savename = ['DSCDL_BID_AN_ADPU_AP_CCNoise/DSCDL_BID_AN_ADPU_AP_CCNoise_' type '_' num2str(nup) '.mat'];
%     save(savename, 'mPSNR', 'mSSIM', 'PSNR', 'SSIM', 'CCmPSNR', 'CCmSSIM', 'CCPSNR', 'CCSSIM');
end