function [im_out, PN] = bscdl_BCGD_ADPU_BID(IMin_y,model,Dict,par,param)
[h,w,ch] = size(IMin_y);
% Used for calculate CSNR
par.height = h;
par.width = w;
par.te_num = ch;
% Initial
im_out = IMin_y;
PN = cell(1,par.cls_num);
for t = 1 : par.nInnerLoop
    if t == 1
        psf = fspecial('gaussian', par.win+2, 2.2);
        YC = data2patch(conv2(im_out, psf, 'same') - im_out, par);
        meanYC = repmat(mean(YC), [par.win^2 1]);
        YC = YC - meanYC;
        AN = zeros(par.K, size(YC, 2));
        AC = zeros(par.K, size(YC, 2));
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
            Alphan = AN(:, idx_cluster);
            Alphac = AC(:, idx_cluster);
        end
        D = [Dn zeros(size(Dc)); Dc zeros(size(Dn)); -par.sqrtmu*Un par.sqrtmu*Uc];
        Y = [Xn - Pn; Xc; zeros(size(Un*Alphan))];
        A12 = [Alphan; Alphac];
        Alpha = zeros(2*size(Dn,2),size(Xn,2));
        for j = 1:size(Xn,2)
            Alpha(:,j) = cgdsq(Y(:,j), D, A12(:,j), param);
            Alphan(:,j) = Alpha(1:size(Alpha,1)/2,j);
            Alphac(:,j) = Alpha(size(Alpha,1)/2+1:end,j);
        end
        Pn = Xn - Dn * Alphan; % 08/17/2016;
        PN{i} = Pn;
        clear Y D A;
        %% Reconstruction
        Xc = Dc * Alphac;
        XC(:, idx_cluster) = Xc;
        AN(:, idx_cluster) = Alphan;
        AC(:, idx_cluster) = Alphac;
    end
    im_out = patch2data(XC+meanX, h, w, 1,par.win, par.step);
end