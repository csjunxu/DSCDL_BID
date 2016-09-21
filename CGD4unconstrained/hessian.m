%This routine evaluates the Hessian diagonal of the convex quadratic function
function y=hess(x)
n=size(x,1);
y=[4*ones(1,n-1)  2];