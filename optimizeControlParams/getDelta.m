function delta = getDelta(thetanew,policyType)
% GETDELTA computes the reward feedback to be used in the policy update for
% training an artificial agent.
%
% INPUTS:
%   thetanew: orientation to evaluate reward
%   policyType: string that specifies whether the policy should be trained
%       on an internal representation of heading, or directly on the visual
%       scene. Can take values 'heading' or 'visual'; default is 'heading'
%

if nargin<2
    policyType = 'heading';
end


theta0 = linspace(0,2*pi,97);
if strcmp(policyType,'heading')==1
    thetaGrid = linspace(0,2*pi,97);
    delta = zeros(size(thetaGrid));
    delta(thetaGrid<=pi) = thetaGrid(thetaGrid<=pi);
    delta(thetaGrid>pi)  = max(thetaGrid(thetaGrid>pi))-thetaGrid(thetaGrid>pi);
    delta = (delta-mean(delta))./10;
elseif strcmp(policyType,'visual')==1
    thetaGrid = linspace(0,pi,49);
    delta = zeros(size(thetaGrid));
    delta(thetaGrid<=pi/2) = thetaGrid(thetaGrid<=pi/2);
    delta(thetaGrid>pi/2)  = max(thetaGrid(thetaGrid>pi/2))-thetaGrid(thetaGrid>pi/2);
    delta = (delta-mean(delta))./5;
    delta = [delta,fliplr(delta(1:end-1))];
else
    error('unrecognized policy type')
end

[~,ii] = min(abs(thetanew-theta0));
delta = delta(ii);
end