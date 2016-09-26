function data= patch2data(patch, h,w,ch, b, s)
data   =  zeros(h,w,ch, 'single');
im_wei   =  zeros(h,w,ch, 'single');
N       =  h-b+1;
M       =  w-b+1;
r     =  1:s:N;
r     =  [r r(end)+1:N];
c     =  1:s:M;
c     =  [c c(end)+1:M];
k          =   0;
for l = 1:ch
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            data(r-1+i,c-1+j,l)  =  data(r-1+i,c-1+j,l) + reshape( patch(k,:)', [length(r) length(c)]);
            im_wei(r-1+i,c-1+j,l)  =  im_wei(r-1+i,c-1+j,l) + 1;
        end
    end
end
data  =  data./(im_wei+eps);