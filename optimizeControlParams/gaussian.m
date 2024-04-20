function f = gaussian(x,mu,sig)
% GAUSSIAN(X,MU,SIG) = 1/sqrt(2*pi*SIG^2)*exp(-(X-MU)^2/(2*SIG^2)) 
% generates a gaussian function with mean MU and standard deviation SIG  
%

f = 1/sqrt(2*pi*sig.^2) .* exp( -(x-mu).^2 / (2*sig.^2) );