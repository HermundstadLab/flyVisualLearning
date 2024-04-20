function f = vonMises(x,mu,k)
% VONMISES(X,MU,K) = exp(K*cos(X-MU))./(2*pi*besseli(0,K)) generates a von 
% Mises profile with mean MU and concentration K
%

f = exp(k*cos(x-mu))./(2*pi*besseli(0,k));
end