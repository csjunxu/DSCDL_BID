% Convex quadratic function 

function y=fnc(x)
% n=size(x,1);    %n is the number of rows of the column vector x.
% y=(sum(x))^2 + norm(x(1:n-1))^2;
y = norm(x,2);
