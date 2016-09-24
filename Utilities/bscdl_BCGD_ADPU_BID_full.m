function [im_out, PN] = bscdl_BCGD_ADPU_BID_full(IMin_y,model,Dict,par,param)
[h,w,ch] = size(IMin_y);
% Used for calculate CSNR
par.height = h;
par.width = w;
par.te_num = ch;

XN = data2patch(IMin_y,  par);
% Initial
im_out = IMin_y;
PN = cell(1,par.cls_num);
AN = zeros(par.K, size(XN, 2));
AC = zeros(par.K, size(XN, 2));
for t = 1 : par.nInnerLoop
    if t == 1
        psf = fspecial('gaussian', par.win+2, 2.2);
        YC = data2patch(conv2(im_out, psf, 'same') - im_out, par);
        meanYC = repmat(mean(YC), [par.win^2 1]);
        YC = YC - meanYC;
        %% GMM: full posterior calculation
        PYZ = zeros(model.nmodels,size(YC,2));
        for i = 1:model.nmodels
            sigma = model.covs(:,:,i);
            [R,~] = chol(sigma);
            Q = R'\YC;
            PYZ(i,:)  = - sum(log(diag(R))) - dot(Q,Q,1)/2;
        end
        %% find the most likely component for each patch group
        [~,cls_idx] = max(PYZ);
    end
    XC = data2patch(im_out,  par);
    XN = data2patch(IMin_y,  par); % one time is ok
    meanX = repmat(mean(XC), [par.win^2 1]);
    XN = XN - meanX;
    XC = XC - meanX;
    for i = 1 : par.cls_num
        idx_cluster   = find(cls_idx == i);
        Xc    = double(XC(:, idx_cluster));
        Xn    = double(XN(:, idx_cluster));
        Dc    = Dict.DC{i};
        Dn    = Dict.DN{i};
        Uc    = Dict.UC{i};
        Un    = Dict.UN{i};
        if (t == 1)
            Pn    = zeros(size(Xn));
            Alphan = mexLasso(Xn - Pn, Dn, param);
            Alphac = Uc \ Un * Alphan;
            Xc = Dc * Alphac; % Xc->Xn; 07/07/2016;  Xn->Xc; 07/22/2016;
        else
            Pn = PN{i};
            Alphac = AC(:, idx_cluster);
        end
        D = [Dn; par.sqrtmu * Un]; % Wn ->Un 07/22/2016;
        Y = [Xn - Pn; par.sqrtmu * Uc * full(Alphac)];
        %  Pn = (Xn - Dn * Alphan) / (1 + par.nu); % 08/17/2016;
        %  Pn can't be put here since there is no Alphan when t > 1.
        Alphan = mexLasso(Y, D,param);
        Pn = Xn - Dn * Alphan; % 08/17/2016;
        PN{i} = Pn;
        clear Y D;
        %% CVPR2012 SCDL case
        D = [Dc; par.sqrtmu * Uc];
        Y = [Xc; par.sqrtmu * Un * full(Alphan)];
        Alphac = full(mexLasso(Y, D,param));
        clear Y D;
        %             %% ICCV2013 MML case
        %             Alphac = Uc \ Un * Alphan;
        %% Reconstruction
        Xc = Dc * Alphac;
        XC(:, idx_cluster) = Xc;
        AN(:, idx_cluster) = Alphan;
        AC(:, idx_cluster) = Alphac;
    end
    im_out = patch2data(XC+meanX, h, w, 1,par.win, par.step);
end