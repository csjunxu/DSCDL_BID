function [NPGset, CPGset,NPG0set, CPG0set] = sample_PG_couple(Nim, Cim, par, patch_num, R_thresh)
if size(Nim, 3) == 3,
    Nim = rgb2ycbcr(Nim);
else
    disp('grayscale image sampled!');
end
[nrow, ncol] = size(Nim);
par.maxr         =  nrow - par.patch_size + 1;
par.maxc         =  ncol - par.patch_size + 1;
r         =  1:par.maxr;
par.r         =  [r r(end)+1:par.maxr];
c         =  1:par.maxc;
par.c         =  [c c(end)+1:par.maxc];
% Index image patch
Index     =   (1:par.maxr*par.maxc);
par.Index    =   reshape(Index, par.maxr, par.maxc);

k    =  0;
CPatches = zeros(par.patch_size^2,length(r)*length(c), 'double');
NPatches = zeros(par.patch_size^2,length(r)*length(c), 'double');
% first j then i : column wise;  first i then j : row wise;
% since N = Npatch(:) is column wise, we employ first j then i.
for j  = 1:par.patch_size
    for i  = 1:par.patch_size
        k         =  k+1;
        Nblk    =  Nim(r-1+i,c-1+j);
        NPatches(k,:) =  Nblk(:)';
        Cblk    =  Cim(r-1+i,c-1+j);
        CPatches(k,:) =  Cblk(:)';
    end
end

x = randperm(nrow-2*par.patch_size-1) + par.patch_size;
y = randperm(ncol-2*par.patch_size-1) + par.patch_size;

[X,Y] = meshgrid(x,y);

xrow = X(:);
ycol = Y(:);

Nim = double(Nim);
Cim = double(Cim);

NPGset = [];
CPGset = [];
NPG0set = [];
CPG0set = [];

idx=1;ii=1;n=length(xrow);
while (idx < patch_num) && (ii<=n),
    row = xrow(ii);
    col = ycol(ii);
    Npatch = Nim(row:row+par.patch_size-1,col:col+par.patch_size-1);
    Cpatch = Cim(row:row+par.patch_size-1,col:col+par.patch_size-1);
    % get PG from P
    NPG = Get_PGfP(NPatches,row,col,par);
    CPG = Get_PGfP(CPatches,row,col,par);
    % check if it is a stochastic patch
    C = Cpatch(:) - mean(Cpatch(:));
    Cnorm=sqrt(sum(C.^2));
    C_normalised=reshape(C/Cnorm,par.patch_size,par.patch_size);
    % pick out small variance patch
    if max(var(CPG)) <= 0.001
        NPG0set = [NPG0set NPG];
        CPG0set = [CPG0set CPG];
        idx=idx+1;
    else
        % eliminate stochastic patch
        if dominant_measure(C_normalised)>R_thresh
            %if dominant_measure_G(Lpatch1,Lpatch2)>R_thresh
            NPGset = [NPGset NPG];
            CPGset = [CPGset CPG];
            idx=idx+1;
        end
    end
    
    ii=ii+1;
end

fprintf('sampled %d patches.\r\n',patch_num);
end

function R = dominant_measure(p)
% calculate the dominant measure
% ref paper: Eigenvalues and condition numbers of random matries, 1988
% p = size n x n patch

hf1 = [-1,0,1];
vf1 = [-1,0,1]';
Gx = conv2(p, hf1,'same');
Gy = conv2(p, vf1,'same');

G=[Gx(:),Gy(:)];
[U, S, V]=svd(G);

R=(S(1,1)-S(2,2))/(S(1,1)+S(2,2));

end
