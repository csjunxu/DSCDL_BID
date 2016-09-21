function im_pout= patch2data(Y, h,w,ch, b, s)
im_pout   =  zeros(h,w,ch, 'single');
im_wei   =  zeros(h,w,ch, 'single');
k          =   0;
N       =  h-b+1;
M       =  w-b+1;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
N       =  length(r);
M       =  length(c);
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        im_pout(r-1+i,c-1+j,:)  =  im_pout(r-1+i,c-1+j,:) + reshape( Y(k,:)', [N M ch]);
        im_wei(r-1+i,c-1+j,:)  =  im_wei(r-1+i,c-1+j,:) + 1;       
    end
end
im_pout  =  im_pout./(im_wei+eps);