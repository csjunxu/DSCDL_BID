clear;
warning off;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');

Original_image_dir = '../../ECCV2016/grayimages/';
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

load Data/params.mat par param;
load Data/Dict_DSCDL_backup_flexible_20160724T132738.mat Dict;
load Data/EMGM_8x8_100_knnNI2BS500Train_20160722T082406.mat;
par.nInnerLoop = 1;
par.cls_num = 100;
type = 'all';
% 'random';
% 'middle';
% 'all';

% Gaussian Noise Removal
sigma = 10;
for sqrtmu = [0.1:0.05:0.5]
    param.lambda = 0.0001;
    par.sqrtmu = sqrtmu;
    PSNR = [];
    SSIM = [];
    for i = 1:im_num
        IMin0=im2double(imread(fullfile(Original_image_dir, im_dir(i).name)));
        S = regexp(im_dir(i).name, '\.', 'split');
        IMname = S{1};
        randn('seed',0);
        IMin = IMin0 + sigma/255*randn(size(IMin0));
        [h,w,ch] = size(IMin);
        hh = [0:1000:h,h];
        ww = [0:1000:w,w];
        num_part = 0;
        if strcmp(type, 'all')
            listh = 1 : length(hh)-1;
            listw = 1 : length(ww)-1;
            IMout = zeros(h,w,ch);
        elseif strcmp(type, 'random')
            listh = randi([1, length(hh)-1],1,1);
            listw = randi([1, length(ww)-1],1,1);
            IMout = zeros(min(h,1000),min(w,1000),ch);
        elseif strcmp(type, 'middle')
            listh = floor(median(1:length(hh)-1));
            listw = floor(median(1:length(ww)-1));
            IMout = zeros(min(h,1000),min(w,1000),ch);
        end
        %%
        for nh = listh
            for nw = listw
                num_part = num_part + 1;
                IMin_part = IMin(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:);
                IMin0_part = IMin0(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:);
                % color or gray image
                if ch==1
                    IMin_part_y = IMin_part;
                    IMin0_part_y = IMin0_part;
                else
                    % change color space, work on illuminance only
                    IMin_part_ycbcr = rgb2ycbcr(IMin_part);
                    IMin_part_y = IMin_part_ycbcr(:, :, 1);
                    IMin_part_cb = IMin_part_ycbcr(:, :, 2);
                    IMin_part_cr = IMin_part_ycbcr(:, :, 3);
                    IMin0_part_ycbcr = rgb2ycbcr(IMin0_part);
                    IMin0_part_y = IMin0_part_ycbcr(:, :, 1);
                    IMin0_part_cb = IMin0_part_ycbcr(:, :, 2);
                    IMin0_part_cr = IMin0_part_ycbcr(:, :, 3);
                end
                fprintf('%s:\n The initial PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name,csnr( IMin_part_y*255, IMin0_part_y*255, 0, 0 ),cal_ssim( IMin_part_y*255, IMin0_part_y*255, 0, 0 ));
                %%
                par.nOuterLoop = 1;
                Continue = true;
                while Continue
                    fprintf('OutIter: %d \n', par.nOuterLoop);
                    IMout_part_y = bscdl_GD_full(IMin_part_y,IMin0_part_y,model,Dict,par,param);
                    % Noise Level Estimation
                    nSig = NoiseLevel(IMout_part_y*255);
                    fprintf('The noise level is %2.4f.\n',nSig);
                    if nSig < 1 || par.nOuterLoop >= 3
                        Continue = false;
                    else
                        par.nOuterLoop = par.nOuterLoop + 1;
                        IMin_part_y = IMout_part_y;
                    end
                end
                    fprintf('The nOurerLoop of the %d/%d part is : %d \n', num_part, (length(hh)-1)*(length(ww)-1), par.nOuterLoop);
                    if ch==1
                        IMout_part = IMout_part_y;
                    else
                        IMout_part_ycbcr = zeros(size(IMin_part));
                        IMout_part_ycbcr(:, :, 1) = IMout_part_y;
                        IMout_part_ycbcr(:, :, 2) = IMin_part_cb;
                        IMout_part_ycbcr(:, :, 3) = IMin_part_cr;
                        IMout_part = ycbcr2rgb(IMout_part_ycbcr);
                    end
                    %% output
                    %                         fprintf('This is the %d/%d part of the image %s.%s!\n',num_part,(length(hh)-1)*(length(ww)-1),IMname,S{2});
                    IMout(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:) = IMout_part;
                end
        end
        %% output
        PSNR = [PSNR csnr( IMout*255,IMin0*255, 0, 0 )];
        SSIM = [SSIM cal_ssim( IMout*255, IMin0*255, 0, 0 )];
        fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name,csnr( IMout*255, IMin0*255, 0, 0 ),cal_ssim( IMout*255, IMin0*255, 0, 0 ));
        imwrite(IMout, ['./DSCDL_GD/DSCDL_GD_AN_' num2str(sigma) '_' IMname '_' num2str(param.lambda) '_' num2str(par.sqrtmu) '.png']);
    end
    %% save output
    mPSNR = mean(PSNR);
    mSSIM = mean(SSIM);
    result = sprintf('./DSCDL_GD/DSCDL_GD_AN_%d_%2.4_%2.2f.mat',sigma,param.lambda,par.sqrtmu);
    save(result,'nSig','PSNR','SSIM','mPSNR','mSSIM');
end