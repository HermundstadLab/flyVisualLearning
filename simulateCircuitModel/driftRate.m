function nu = driftRate(theta,heading,w,nx)
% DRIFTRATE returns the output of PFL2 neurons as a function of heading 
% relative to a goal. 
%
% INPUTS:
%   theta: array of all possible heading values 
%   heading: current compass heading
%   w:       weight vector that defines goal heading
%   nx:      number of unique heading values (should be numel(theta))
%
% OUTPUT:
%   nu: drift rate
%

[rO,B] = outputNeuron(theta,heading,w,nx,'center');
nu = rO+B;