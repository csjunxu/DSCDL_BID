function [PG, PGmean] = sample_PatchGroups(im, patch_num, par)
if size(im, 3) == 3,
    im = rgb2ycbcr(im);
end

if par.Patch_Channel == 3
    disp('RGB image sampled!');
else
    im = im(:,:,1);
    disp('Grayscale image sampled!');
end

[nrow, ncol, nch] = size(im);
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
Patches = zeros(par.patch_size^2,length(r)*length(c), 'double');
% first j then i : column wise;  first i then j : row wise;
% since N = Npatch(:) is column wise, we employ first j then i.
for l = 1:nch
    for j  = 1:par.patch_size
        for i  = 1:par.patch_size
            k         =  k+1;
            blk    =  im(r-1+i,c-1+j,l);
            Patches(k,:) =  blk(:)';
        end
    end
end
x = randperm(nrow-2*par.patch_size-1) + par.patch_size;
y = randperm(ncol-2*par.patch_size-1) + par.patch_size;

[X,Y] = meshgrid(x,y);

xrow = X(:);
ycol = Y(:);

im = double(im);
im = double(im);

PG = [];
PGmean = [];
idx=1;ii=1;n=length(xrow);
while (idx < patch_num) && (ii<=n),
    row = xrow(ii);
    col = ycol(ii);
    patch = im(row:row+par.patch_size-1,col:col+par.patch_size-1,:);
    % get PG from P
    [PatchGroup, mean] = Get_PGfP(Patches,row,col,par);
    PG = [PG PatchGroup];
    PGmean = [PGmean mean];
    idx=idx+1;
    ii=ii+1;
end
fprintf('sampled %d patches.\r\n',patch_num);
end
