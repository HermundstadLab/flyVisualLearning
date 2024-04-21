function [headLocs,goalLocs,headGains,goalGains,PIscores,distDanger,angdiff,headingVisits,headingJumps,gain] = sweepPIscores(p,params)
% SWEEPPISCORES runs many simulations of an agent (model fly) that must form
% a mapping of the visual environment while learning to avoid danger zones.
% Simulations vary in the initial strength of the goal heading. Results are
% saved in directory sims/.
% 
% INPUTS:
%   p:      structure containing simulation parameters
%   params: cell array that specifies simulation setting
%       params{1}: 'asym' or 'sym' (whether current scene is asymmetric or symmetric)
%       params{2}: 'train' or 'probe' (whether or not model flies experience punishment)
%       params{3}: 'uniform' or 'random' (whether initial goal location is uniform or random across model flies)
%       params{4}: 'sym_init' or 'asym_init' (whether visual map is initialized in a symmetric or asymmetric scene)
%
% OUTPUTS:
%   headLocs:       location of most stable compass heading over time
%   goalLocs:       location of goal heading over time
%   headGains:      strength of most stable compass heading over time
%   goalGains:      strength of goal heading over time
%   PIscores:       PI score over time
%   distDanger:     angular distance between goal heading and center of safe 
%                       zone over time 
%   angDiff:        angular distance between most stable compass heading and
%                       goal heading
%   headingVisits:  histogram of visits to different compass headings
%   headingJumps:   histogram of headings at which compass heading jumped
%   gain:           array containing initial strengths of goal heading
%

if nargin<2
    params =  {'sym','train','uniform','asym_init'};
end

%initial conditions 
gain = linspace(0,.5,6);
gain = gain(2:end);
niter = 200;
centerDangerLocs = [1,17];

if strcmp(params{2},'train')==1 
    if strcmp(params{3},'uniform')==1
        loc = centerDangerLocs(1).*ones(1,niter);    % center goal heading within danger zone
    end
elseif strcmp(params{2},'probe')==1
    if strcmp(params{3},'uniform')==1
        loc = centerDangerLocs(1).*ones(1,niter);    % center goal heading within danger zone
    end
else
    error('unrecognized simulation parameters')
end

if strcmp(params{3},'random')==1
    locs  = randi(p.N,[1,5*niter]);
    gains = rand(1,5*niter).*0.5;
else
    locs  = repmat(loc,[1,numel(gain)]);
    gains = [];
    for i=1:numel(gain)
        gains = [gains,repmat(gain(i),[1,niter])];
    end
end


headLocs   = [];
goalLocs   = [];
headGains  = [];
goalGains  = [];
PIscores   = [];
distDanger = [];
angdiff    = [];
gain       = [];
headingVisits = [];
headingJumps  = [];

parfor i=1:numel(gains)
    [headLocsTmp,goalLocsTmp,headGainsTmp,goalGainsTmp,...
        PIscoresTmp,distDangerTmp,angdiffTmp,headingVisitsTmp,headingJumpsTmp,~]...
        = runLearningSims(p,gains(i),locs(i),params);
    
    headLocs   = [headLocs;  headLocsTmp  ];
    goalLocs   = [goalLocs;  goalLocsTmp  ];
    headGains  = [headGains; headGainsTmp ];
    goalGains  = [goalGains; goalGainsTmp ];
    PIscores   = [PIscores;  PIscoresTmp  ];
    distDanger = [distDanger;distDangerTmp];
    angdiff    = [angdiff;   angdiffTmp   ];
    gain          = [gain; repmat(gains(i),[niter,1])];
    headingVisits = [headingVisits; headingVisitsTmp ];
    headingJumps  = [headingJumps;  headingJumpsTmp  ];
    
end


save(['sims/PIsweep_',params{1},'_',params{2},'_',params{3},'_',params{4},'.mat'],'headLocs','goalLocs','headGains',...
    'goalGains','PIscores','distDanger','angdiff','headingVisits','headingJumps','gain')
end


