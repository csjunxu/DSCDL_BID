function im_out = bscdl_BID_full(IMin_y,model,Dict,par,param)
[h, w, ch] = size(IMin_y);
% Initial
im_out = IMin_y;
XN = data2patch(IMin_y,  par);
AN = zeros(par.K, size(XN, 2));
AC = zeros(par.K, size(XN, 2));
for t = 1 : par.nInnerLoop
    if t == 1
        psf = fspecial('gaussian', par.win+2, 2.2);
        YH = data2patch(conv2(im_out, psf, 'same') - im_out, par);
        meanY = repmat(mean(YH), [par.win^2 1]);
        YH = YH - meanY;
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
    XC = data2patch(im_out,  par);
    XN = data2patch(IMin_y,  par);
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
        %% Reconstruction
        Xc = Dc * Alphac;
        XC(:, idx_cluster) = Xc;
        AN(:, idx_cluster) = Alphan;
        AC(:, idx_cluster) = Alphac;
    end
    im_out = patch2data(XC+meanX, h, w, 1,par.win, par.step);
    [N, ~]       =   Compute_NLM_Matrix( im_out , 5, par);
    NTN          =   N'*N*0.05;
    im_f = sparse(double(im_out(:)));
    for i = 1 : fix(60 / t.^2)      
        im_f = im_f  - NTN*im_f;
    end
    im_out = reshape(full(im_f), h, w);
end