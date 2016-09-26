function patch = data2patch(data, par)
[h, w, ch]   =   size(data);
b          =   par.win;
s          =   par.step;
N       =  h-b+1;
M       =  w-b+1;
r     =  1:s:N;
r     =  [r r(end)+1:N];
c     =  1:s:M;
c     =  [c c(end)+1:M];
patch      =  zeros(b*b*ch,length(r)*length(c), 'single');
k          =   0;
for l = 1:ch
    for i  = 1:b
        for j  = 1:b
            k        =  k+1;
            blk  =  data(r-1+i,c-1+j,l);
            patch(k,:)  =  blk(:)';
        end
    end
end