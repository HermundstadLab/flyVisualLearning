function x = nonlin(x)
% NONLIN(X) returns the heaviside function of X
%

x(x<0)=0;
x(x>0)=1;
end