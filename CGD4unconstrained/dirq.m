%  dirq.m   (May 20, 2006 (revised Oct 10, 2006))
%
%This routine computes the direction d whose ith component equals
%that of q_H=(-q_H(x;1),...,-q_H(x;n))^T if the latter is greater than 
%upsilon*||q_H||_infty in magnitude, otherwise equals zero.
%Here H is a diagonal positive definite matrix.
%It also outputs maxRr, which is ||Hd_{H}(x)||_{\infty} and will be used
%for checking the termination.

function [maxRr,d] = dirq(c,x,g,h,ups)

R=-median([ x' ; (g'+c)./h ; (g'-c)./h ]);     % R=d_{H}(x)
absR=abs(R);
maxR=max(absR);
hR=h.*R;
Q=-g.*R'-.5*R'.*hR'-c*abs(x+R')+c*abs(x);      % Q=q_H
%sumQ=sum(Q);                                   % sumQ=||q_H(x)||_1
maxQ=max(Q);                                   % ||q_H||_infty

maxRr=max(abs(hR));                            % maxRr=||Hd_{H}(x)||_{\infty} 

%indx=find(Q'<ups*sumQ);
indx=find(Q'<ups*maxQ);
d=R';
d(indx)=0;	                               % set d(i)=0 if Q(i)<ups*maxQ