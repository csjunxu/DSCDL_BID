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
load Data/params_flexible.mat par param;
load Data/DSCDL_BID_Dict_ADPU_backup_nup0.5_20160819T042830.mat Dict;
load Data/EMGM_8x8_100_knnNI2BS500Train_20160722T082406.mat;
par.cls_num = 100;
par.nInnerLoop = 3;
%% the whole image or part
type = 'middle';
% 'random';
% 'middle';
% 'all';
for i = 1 : im_num
    IMin=im2double(imread(fullfile(Original_image_dir, im_dir(i).name)));
    S = regexp(im_dir(i).name, '\.', 'split');
    IMname = S{1};
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
            imwrite(IMin_part, ['./TestedImages/' IMname '_' type '.png']);

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
%             while Continue
%                 fprintf('Iter: %d \n', nOuterLoop);
%                 IMout_part_y = bscdl_BID_full(IMin_part_y,model,Dict,par,param);
%                 % Noise Level Estimation
%                 nSig = NoiseLevel(IMout_part_y*255);
%                 fprintf('The noise level is %2.4f.\n',nSig);
%                 if nSig < 0.005 || nOuterLoop >= 10
%                     Continue = false;
%                 else
%                     nOuterLoop = nOuterLoop + 1;
%                     IMin_part_y = IMout_part_y;
%                 end
%             end
%             if ch==1
%                 IMout_part = IMout_part_y;
%             else
%                 IMout_part_ycbcr = zeros(size(IMin_part));
%                 IMout_part_ycbcr(:, :, 1) = IMout_part_y;
%                 IMout_part_ycbcr(:, :, 2) = IMin_part_cb;
%                 IMout_part_ycbcr(:, :, 3) = IMin_part_cr;
%                 IMout_part = ycbcr2rgb(IMout_part_ycbcr);
%             end
%             fprintf('The nOurerLoop of the %d/%d part is : %d \n', num_part, (length(hh)-1)*(length(ww)-1), nOuterLoop);
%             if strcmp(type, 'all')
%                 IMout(hh(nh)+1:hh(nh+1),ww(nw)+1:ww(nw+1),:) = IMout_part;
%             elseif strcmp(type, 'random') || strcmp(type, 'middle')
%                 IMout = IMout_part;
%             end
        end
    end
    %% output
%     imwrite(IMout, ['./DSCDL_BID_AN_ADPU/ADPU_nup0.5_DSCDL_BID_AN_' type '_' IMname '.png']);
%     fprintf('The Loops of these parts are  %s. \n', num2str(nOuterLoop));
%     Loops = nOuterLoop*par.nInnerLoop;
%     save(['./DSCDL_BID_AN_ADPU/ADPU_nup0.5_DSCDL_BID_AN_' IMname '_Loops.mat'],'Loops');
end