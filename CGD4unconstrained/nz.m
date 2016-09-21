% May 20, 2006 (revised Oct 04, 2006)
%  nz.m 
%This routine computes the number of components of x whose absolute value
%is greater than and equal to t. 

function z  = nz(x,n)

t=1e-15;
xx=zeros(n,1);
xx(find(abs(x)>=t))=1;
z=sum(xx);

%xx=x;
%xx(find(t>abs(x)))=0;
%xx(find(t<=abs(x)))=1;
%z=xx'*ones(n,1);