clear;
warning off;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
addpath('GatingNetwork');

Original_image_dir = './TestingImages/';
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

cls_num = 100;
for rho = 0.05
    par.rho = rho;
    params = 'Data/params__flexible.mat';
    load(params,'par','param');
    Dict_SR_backup = 'Data/Dict_SR_backup_flexible.mat';
    load(Dict_SR_backup,'Dict');
    %% load pre-trained models (see demoTrain.m for training models)
    model = getfield(load('./GatingNetwork/tmpGMM.mat'),'model');
    model.net = getfield(load('./GatingNetwork/tmpGatingNet.mat'),'net');
    par.nInnerLoop = 5;
    
    for i = 1:im_num
        par.nOuterLoop = 1;
        IMin=im2double(imread(fullfile(Original_image_dir, im_dir(i).name)));
        S = regexp(im_dir(i).name, '\.', 'split');
        IMname = S{1};
        fprintf('The Input Image is %s.\n',IMname);
        % color or gray image
        [h,w,ch] = size(IMin);
        if h >= 600
            IMin = IMin(ceil(h/2)-300+1:ceil(h/2)+300,:,:);
        end
        if w >= 800
            IMin = IMin(:,ceil(w/2)-400+1:ceil(w/2)+400,:);
        end
        [h,w,ch] = size(IMin);
        imwrite(IMin, ['Input_' IMname '.png']);
        if ch==1
            IMin_y = IMin;
        else
            % change color space, work on illuminance only
            IMin_ycbcr = rgb2ycbcr(IMin);
            IMin_y = IMin_ycbcr(:, :, 1);
            IMin_cb = IMin_ycbcr(:, :, 2);
            IMin_cr = IMin_ycbcr(:, :, 3);
        end
        Continue = true;
        while Continue
            fprintf('Iter: %d \n', par.nOuterLoop);
            nSig = NoiseLevel(IMin_y*255);
            fprintf('The noise level before denoising is %2.2f.\n',nSig);
            tic
            [Iout_y] = scdl_interp(IMin_y,model,Dict,par,param);
            toc
            % Noise Level Estimation
            nSig = NoiseLevel(Iout_y*255);
            fprintf('The noise level after denoising is %2.2f.\n',nSig);
            Iout_y(Iout_y>1)=1;
            Iout_y(Iout_y<0)=0;
            if ch==1
                Iout = Iout_y;
            else
                Iout_ycbcr = zeros([h,w,ch]);
                Iout_ycbcr(:, :, 1) = Iout_y;
                Iout_ycbcr(:, :, 2) = IMin_cb;
                Iout_ycbcr(:, :, 3) = IMin_cr;
                Iout = ycbcr2rgb(Iout_ycbcr);
            end
            %% output
            imwrite(Iout, ['./Denoised_flexible_NE_GN/Denoised_flexible_' IMname '_NE_GN_Ori_' num2str(par.rho) '_nInnerLoop_' num2str(par.nInnerLoop) '_nOuterLoop_' num2str(par.nOuterLoop) '.png']);
            if nSig < 0.1 || par.nOuterLoop >= 10
                Continue = false;
            else
                par.nOuterLoop = par.nOuterLoop + 1;
                IMin_y = Iout_y;
            end
        end
    end
end