function im_out = DSCDL_PG_BID_20161001(IMin_y,IM_GT_y,model,Dict,par,param)
%% Initialization
im_out = IMin_y;
for t = 1 : par.nInnerLoop
    if t == 1
        psf = fspecial('gaussian', par.win+2, 2.2);
        [nDCnlYH,~,~,par] = Image2PGs( conv2(im_out, psf, 'same') - im_out, par);
        AN = zeros(par.K, size(nDCnlYH, 2));
        AC = zeros(par.K, size(nDCnlYH, 2));
        %% GMM: full posterior calculation
        nPG = size(nDCnlYH,2)/par.nlsp; % number of PGs
        PYZ = zeros(model.nmodels,nPG);
        for i = 1:model.nmodels
            sigma = model.covs(:,:,i);
            [R,~] = chol(sigma);
            Q = R'\nDCnlYH;
            TempPYZ = - sum(log(diag(R))) - dot(Q,Q,1)/2;
            TempPYZ = reshape(TempPYZ,[par.nlsp nPG]);
            PYZ(i,:) = sum(TempPYZ);
        end
        %% find the most likely component for each patch group
        [~,cls_idx] = max(PYZ);
        cls_idx=repmat(cls_idx, [par.nlsp 1]);
        cls_idx = cls_idx(:);
        [idx,  s_idx] = sort(cls_idx);
        idx2 = idx(1:end-1) - idx(2:end);
        seq = find(idx2);
        seg = [0; seq; length(cls_idx)];
    end
    %%  Image to PGs
    [nDCnlXC,blk_arrXC,DCXC,par] = Image2PGs( im_out, par);
    [nDCnlXN,~,~,par] = Image2PGs( IMin_y, par);
    X_hat = zeros(par.ps^2,par.maxr*par.maxc,'double');
    W = zeros(par.ps^2,par.maxr*par.maxc,'double');
    for   i  = 1:length(seg)-1
        idx          =   s_idx(seg(i)+1:seg(i+1));
        cls       =   cls_idx(idx(1));
        Xc    = nDCnlXC(:, idx);
        Xn    = nDCnlXN(:, idx);
        Dc    = Dict.DC{cls};
        Dn    = Dict.DN{cls};
        Uc    = Dict.UC{cls};
        Un    = Dict.UN{cls};
        if (t == 1)
            Alphan = mexLasso(Xn, Dn, param);
            Alphac = Uc \ Un * Alphan;
            Xc = Dc * Alphac;
        else
            Alphac = AC(:, idx_cluster);
        end
        %% Transformation
        D = [Dn; par.sqrtmu * Un]; 
        Y = [Xn; par.sqrtmu * Uc * full(Alphac)];
        Alphan = mexLasso(Y, D,param);
        clear Y D;
        D = [Dc; par.sqrtmu * Uc];
        Y = [Xc; par.sqrtmu * Un * full(Alphan)];
        Alphac = full(mexLasso(Y, D,param));
        clear Y D;
        %% Reconstruction
        Xc = Dc * Alphac;
        nDCnlXC(:, idx) = Xc;
        AN(:, idx) = Alphan;
        AC(:, idx) = Alphac;
        X_hat(:,blk_arrXC(idx)) = X_hat(:,blk_arrXC(idx)) + nDCnlXC(:, idx) + DCXC(:,idx);
        W(:,blk_arrXC(idx))     = bsxfun(@plus,W(:,blk_arrXC(idx)),ones(par.ps^2,1));
    end
    %% PGs to Image
    im_out = PGs2Image(X_hat,W,par);
    fprintf('nInnerLoop: %d, PSNR = %2.4f, SSIM = %2.4f. \n', t, csnr( im_out*255, IM_GT_y*255, 0, 0 ),cal_ssim( im_out*255, IM_GT_y*255, 0, 0 ));
end