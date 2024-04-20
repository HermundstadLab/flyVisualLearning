function f = sigmoid(x,k,fmax,offset)
% SIGMOID(X,K,FMAX,OFFSET) = FMAX/(1+exp(-K*X))-OFFSET generates a 
% sigmoidal function with slope K, scale FMAX , and vertical shift OFFSET. 
%

f = fmax./(1+exp(-k*x)) - offset;
end