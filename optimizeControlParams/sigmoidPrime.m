function f = sigmoidPrime(x,k,fmax)
% SIGMOIDPRIME(X,K,FMAX) = FMAX*K*exp(-K*X))/((1+exp(-K*X))^2 evaluates the 
% derivative of the sigmoid function (ith slope K, scale FMAX) at X. 
%

f = (fmax*k*exp(-k*x))./((1+exp(-k*x)).^2);
end