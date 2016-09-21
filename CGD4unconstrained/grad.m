%This routine evaluates the gradient of the convex quadratic function
function y=grad(x)
n=size(x,1);
y=2*( ones(n,1)*sum(x) + [x(1:n-1);0] );

