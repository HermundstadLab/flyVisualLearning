function pR = turnBias(theta,heading,w,nx)
% TURNBIAS returns the output of PFL3 neurons as a function of heading 
% relative to a goal. 
%
% INPUTS:
%   theta: array of all possible heading values 
%   heading: current compass heading
%   w:       weight vector that defines goal heading
%   nx:      number of unique heading values (should be numel(theta))
%
% OUTPUT:
%   pR: probability of rightward saccade

[rR,BR] = outputNeuron(theta,heading,w,nx,'right');
[rL,BL] = outputNeuron(theta,heading,w,nx,'left');

pR = .5*( 1 + (rR+BR) - (rL+BL) );