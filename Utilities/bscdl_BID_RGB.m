function IMout = bscdl_BID_RGB(IMin,model,Dict,par,param)
[h, w, ch] = size(IMin);
% Initial
IMout = IMin;
for t = 1 : par.nInnerLoop
    if t == 1
        psf = fspecial('gaussian', par.win+2, 2.2);
        YH = data2patch(convn(IMout, psf, 'same') - IMout, par);
        meanY = repmat(mean(YH), [par.L 1]);
        YH = YH - meanY;
        AN = zeros(par.K, size(YH, 2));
        AC = zeros(par.K, size(YH, 2));
        %% GMM: full posterior calculation
        PYZ = zeros(model.nmodels,size(YH,2));
        for i = 1:model.nmodels
            sigma = model.covs(:,:,i);
            [R,~] = chol(sigma);
            Q = R'\YH;
            PYZ(i,:)  = - sum(log(diag(R))) - dot(Q,Q,1)/2;
        end
        %% find the most likely component for each patch group
        [~,cls_idx] = max(PYZ);
    end
    XC = data2patch(IMout,  par);
    XN = data2patch(IMin,  par);
    meanX = repmat(mean(XC), [par.L 1]);
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
            Alphan = mexLasso(Xn, Dn, param);
            Alphac = Uc \ Un * Alphan;
            Xc = Dc * Alphac;
        else
            Alphac = AC(:, idx_cluster);
        end
        D = [Dn; par.sqrtmu * Un]; % Wn ->Un 07/22/2016;
        Y = [Xn; par.sqrtmu * Uc * full(Alphac)];
        Alphan = mexLasso(Y, D,param);
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
    IMout = patch2data(XC+meanX, h, w, ch, par.win, par.step);
end