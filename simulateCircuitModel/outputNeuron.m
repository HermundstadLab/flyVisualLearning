function [rO,B] = outputNeuron(theta,heading,w,nx,side)
% OUTPUTNEURON returns the output of PFL 'action' neurons as a function of
% heading relative to a goal. 
%
% INPUTS:
%   theta: array of all possible heading values 
%   heading: current compass heading
%   w:       weight vector that defines goal heading
%   nx:      number of unique heading values (should be numel(theta))
%   side:    string that specifies whether this output corresponds to the
%               right/left PFL3 population, or the PFL2 population. Can
%               take values 'right', 'left', or 'center'
%

if strcmp(side,'right')==1
    ind = heading-nx/4;
    ind = circIndex(ind,nx);
    p   = firingRate(cos(theta-theta(ind)));
    rO  = sum(w.*p)./sum(p.*p);
    B = 0;
elseif strcmp(side,'left')==1
    ind = heading+nx/4;
    ind = circIndex(ind,nx);
    p   = firingRate(cos(theta-theta(ind)));
    rO  = sum(w.*p)./sum(p.*p);
    B = 0;
elseif strcmp(side,'center')==1
    ind = heading-nx/2;
    ind = circIndex(ind,nx);
    p   = firingRate(cos(theta-theta(ind)));
    rO  = sum(w.*p)./sum(p.*p);
    
    vmax = 1.1;
    sumW = sum(w);
    D = sum(p.*p);
    B = vmax-sumW/(D);
    
end
end