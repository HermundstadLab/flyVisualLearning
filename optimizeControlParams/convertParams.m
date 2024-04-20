function [nuF,pR,omegaF,omegaR] = convertParams(omega,fk,params)
% CONVERTPARAMS converts a vector of weights into behavioral control 
% parameters (drift rates and turn biases)
%
% INPUTS:
%   omega: vector of weights output from the learning process
%   fk: a set of basis functions
%   params: input parameters for the sigmoid function that transforms
%       weights into behavior parameters
%
% See also: GETBASISFUNCTIONS
%

nb = (numel(omega))/2;
omegaF = omega(1:nb)';
omegaR = omega(nb+1:2*nb)';

nuF = sigmoid(omegaF*fk,params(1),params(2),params(3));
pR  = sigmoid(omegaR*fk,params(4),params(5),params(6));
end