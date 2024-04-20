function [medDur,locP,scaleP,locD,scaleD] = approxSaccadeParams(sL,sNL)
% APPROXSACCADEPARAMS determines the best-fitting parameters to approximate
% the joint distribution of angular velocities and durations as a single
% distribution of saccade sizes.
%
% INPUTS:
%   sL, sNL:  data structures containining segmented data
%
% OUTPUTS:
%   medDur: median saccade duration
%   locP:   location of best-fitting log-normal distribution (radians)
%   scaleP: scale of best-fitting log-normal distribution (radians)
%   locD:   location of best-fitting log-normal distribution (degrees)
%   scaleD: scale of best-fitting log-normal distribution (degrees)
%

dx = [];
dt = [];

for i=1:numel(sL.fly)
    for j=1:numel(sL.fly(i).trial)
        dx = [dx,abs(sL.fly(i).trial(j).saccades.dx)];
        dt = [dt,(sL.fly(i).trial(j).saccades.dt)];
    end
end
for i=1:numel(sNL.fly)
    for j=1:numel(sNL.fly(i).trial)
        dx = [dx,abs(sNL.fly(i).trial(j).saccades.dx)];
        dt = [dt,(sNL.fly(i).trial(j).saccades.dt)];
    end
end

medDur = median(dt);

%properties of saccade size distribution (computed in pixels):
pd     = fitdist(dx','lognormal');
locP   = pd.mu;
scaleP = pd.sigma;

%converted to degrees:
pd     = fitdist((dx').*360./96,'lognormal');
locD   = pd.mu;
scaleD = pd.sigma;