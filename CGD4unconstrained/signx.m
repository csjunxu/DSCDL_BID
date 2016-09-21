% May 20, 2006 (revised Oct 04, 2006)
%  signx.m 
%
% This routine computes H*d_H(x), t, and the vector s whose 
% ith component is 1 if x_j>=t or -1 if x_j<=-t or 0 otherwise,
% where t = -const/log(min(.001,||R_H(x)||_infty)).

function [s,t,maxRr] = signx(c,x,g,h)

R=-median([ x' ; (g'+c)./h ; (g'-c)./h ]);  % R=d_{H}(x)
absR=abs(R);
maxR=max(absR);
maxRr=max(h.*absR);                         % ||H*d_H(x)||_{\infty}                 
if maxRr==0
  t=0;
else
  t=-.0001/log(min(.1,.01*maxR));
end
s=sign(x);
s(find(t>abs(x)))=0;

