clear;clc;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');
% addpath('SPAMS/release/mkl64');

load Data/params_gray_PG.mat;
task = 'BID';

for j = 1:par.Patch_Channel
    modelname = sprintf('Data/GMM_PG_%d_10_8x8_64_20161003T094301.mat',j);
    eval(['load ' modelname]);
    Dini = cell(par.Patch_Channel,par.cls_num);
    for i = 1 : par.cls_num
        XN_t = double(Xn{j,i});
        XC_t = double(Xc{j,i});
        XN_t = XN_t - repmat(mean(XN_t), [par.win^2 1]);
        XC_t = XC_t - repmat(mean(XC_t), [par.win^2 1]);
        fprintf('DSCDL: Channel: %d, Cluster: %d\n', j, i);
        D = mexTrainDL([XN_t;XC_t], param);
        Dini{j,i} = D;
        save('Data/DSCDL_Dict_PG_3_10_8x8_64_BID_Initial.mat','Dini');
        D = Dini{j,i};
        clear Dini;
        Dn = D(1:par.win * par.win,:);
        Dc = D(par.win * par.win+1:end,:);
        Wn = eye(size(Dn, 2));
        Wc = eye(size(Dc, 2));
        Alphac = mexLasso([XN_t;XC_t], D, param);
        Alphan = Alphac;
        clear D;
        [Alphac, Alphan, XC_t, XN_t, Dc, Dn, Wc, Wn, Uc, Un, Pn, f] = ADPU_Double_Semi_Coupled_DL(Alphac, Alphan, XC_t, XN_t, Dc, Dn, Wc, Wn, par);
        DSCDL.DC{j,i} = Dc;
        DSCDL.DN{j,i} = Dn;
        DSCDL.WC{j,i} = Wc;
        DSCDL.WN{j,i} = Wn;
        DSCDL.UC{j,i} = Uc;
        DSCDL.UN{j,i} = Un;
        DSCDL.PN{j,i} = Pn;
        DSCDL.f{j,i} = f;
        Dict_BID = sprintf('Data/DSCDL_ADPU_Dict_PG_3_10_8x8_64_%s_20161003.mat',task);
        save(Dict_BID,'DSCDL');
    end
end