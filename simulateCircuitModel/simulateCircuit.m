function res = simulateCircuit(p,r)
% SIMULATECIRCUIT is the primary script for running a single-agent
% simulation whose compass heading is shaped by an evolving visual map, and
% whose actions are determined by fixed-form policy tethered to the 
% difference between the compass heading and a flexible goal heading
%
% INPUTS:
%   p:   structure containing simulation parameters
%   r:   structure containing initial conditions
%
% OUTPUTS:
%   res: structure containing simulation results
%

%extract learning rates
alphaW = p.alphaW;
alphaH = p.alphaH;

%sample turn from lognormal distribution
muLN  =  2.56;       % (defined in pixels; muLN = 3.89 in deg)
sigLN = .56;         % (defined in pixels; sigLN = 0.54 in deg)

%assume linear relationship btw turn size and duration
dtscale = 4;         % dt = dtscale*turn (sampled in 100ms time bins)

%partition heading space
nx = p.N;
theta  = linspace(0,2*pi,nx+1);
theta  = theta(1:end-1);
dtheta = theta(2)-theta(1);

pJumpBase = 1;

%-------------------------- for nx = 16 ----------------------------------%
% 2pi  = 16 units
% pi   = 8  units
% pi/2 = 4  units

% distance to goal ranges from -nx/2=-8 to nx/2-1=7
% fixed output curves should be centered on index nx/2+1 = 9
%-------------------------------------------------------------------------%

%------------------------ reward structure -------------------------------%
rD = -1;
rS = +1;
danger     = [rD*ones(1,nx/8),rS*ones(1,nx/4),rD*ones(1,nx/4),rS*ones(1,nx/4),rD*ones(1,nx/8)];
dangerBase = [rD*ones(1,nx/8),rS*ones(1,nx/4),rD*ones(1,nx/4),rS*ones(1,nx/4),rD*ones(1,nx/8)];

%initialize heading
headingLabels = 1:nx;
arenaLabels   = fliplr(headingLabels);
[~,isort] = sort(arenaLabels);

%------------------------ environment types ------------------------------%
% asym
% sym
% symNoJump;

%--------------------------- trial types ---------------------------------%
% frozen (no update to either heading or goal weights)
% probe (update heading weights only)
% train (update heading and goal weights)


if nargin<2
    %goal weights
    w = .3*(-.2*cos(theta-theta(1))-.4*cos(theta-theta(5)).^2+1);

    %heading weights
    hMat = .25*rand(nx)+.5;
    
    %initial headings
    heading = randperm(nx,1);   
    
    trialType = 'train';
    envType   = 'sym';
    
    %initialize training times
    dtTot = 2000;
    mmax  = 20000;
else
    w = r.w;
    hMat = r.hMat;
    heading = r.heading;
    trialType = r.trialType;
    dtTot = r.dtTot;
    mmax  = r.mmax;
    envType = r.envType;
end
arenaHeading = arenaLabels(heading);

%initalize outputs
hMatAll = hMat;
sizes   = nan;
dirs    = nan;
durs    = nan;
jumps   = nan;
rewards = nan;
safeVisits = nan;
arenaHeadings = arenaHeading;
headings      = heading;
bumpProfileIter = [];
wAllIter = w;

visSceneTraj     = [];
arenaHeadingTraj = [];
headingTraj      = [];
indsTraj         = [];
jumpTraj         = [];
offsets          = 0;

% i: indexes discrete actions
% m: indexes time (within actions)
m = 1;
fmax = 1.5; %scales fixation duration
for i = 1:dtTot

    
    if m>mmax
        break
    end
   
    %------------------------- BEGIN ITERATION ---------------------------%
    
    
    %---------------------determine current heading-----------------------%
   
    headings(i)      = heading;
    arenaHeadings(i) = arenaHeading;
      
    %--------------determine punishment at current heading----------------%
    
    if strcmp(trialType,'train')
        R = danger(arenaHeading);
    else
        R = 0;
    end
    
    rewards(i)    = R;
    safeVisits(i) = dangerBase(arenaHeading);
    
    %-------------------determine duration of fixation--------------------%
    nu      = driftRate(theta,heading,w,nx);
    rCout   = fmax./nu;
    durs(i) = rCout;
    durF    = round(rCout*10); %break into 100ms time bins
    
    %----------update goal weights continuously during fixation-----------%
    
    %activity profile at current heading (this will be used to update the 
    %goal (and later heading) weights 
    
    rC = centerHeadingNeuronFB(theta,heading);
    
    for k=1:durF

        if R>0
            dw = relu(rC-w).*nonlin(1-w) - relu(w-rC).*nonlin(w);
        elseif R<0
            dw = -relu(rC-w).*nonlin(w) + relu(w-rC).*nonlin(1-w);
        else
            dw = 0;
        end

        wnew = w+alphaW.*dw;
        w = wnew;

        m=m+1;
        
        visSceneTraj     = [visSceneTraj,isort(arenaHeading)];
        arenaHeadingTraj = [arenaHeadingTraj,arenaHeading];
        headingTraj      = [headingTraj,heading];
        indsTraj         = [indsTraj,1];
        jumpTraj         = [jumpTraj,0];
        
    end
    bumpProfileIter = [bumpProfileIter;rC];
    wAllIter = [wAllIter;w];
    iterNum(i) = m;
    
    
    %-----------------------determine size of saccade---------------------%
    
    turnSize = lognrnd(muLN,sigLN)*(2*pi/96);
    sizes(i) = turnSize;
    
    %--------------------determine direction of saccade-------------------%
    pR = turnBias(theta,heading,w,nx);
    if rand()<pR
        dirR = 1;
        dirL = 0;
    else
        dirL = 1;
        dirR = 0;
    end
    turn = (dirR-dirL)*turnSize;
    turnDur = .3;%dtscale*turnSize;
    dirs(i) = sign(dirR-dirL);
    pRight(i) = pR;
    
    %turns move bump in opposite direction
    if sign(turn)>0
        steps = -1:-1:round(-turn./dtheta);
    else
        steps = 1:1:round(-turn./dtheta);
    end
    ddt = abs(round(-turn./dtheta));
    
    %----------update heading weights continuously during saccade---------%
    headingTemp = heading;
    arenaHeadingTemp = arenaHeading;
    
    %for symmetric environment only:
    %at any given time, two sets of ring neurons are active, and thus the
    %weights at two "candidate arena headings" are updated  
    if strcmp(envType,'sym') || strcmp(envType,'symNoJump')
        candArenaHeadings = [arenaHeading,circIndex(arenaHeading+nx/2,nx)];
    else
        candArenaHeadings =  arenaHeading;
    end
        
    
    for k=1:numel(steps)
        headingTemp      = circIndex(heading     +steps(k),nx);
        arenaHeadingTemp = circIndex(arenaHeading-steps(k),nx);
        
        visSceneTraj     = [visSceneTraj,isort(arenaHeadingTemp)];
        arenaHeadingTraj = [arenaHeadingTraj,arenaHeadingTemp];
        headingTraj      = [headingTraj,headingTemp];
        indsTraj         = [indsTraj,0];
        jumpTraj         = [jumpTraj,0];
        if ~strcmp(trialType,'frozen')
            for mindex=1:numel(candArenaHeadings)

                arenaHeadingUpdate = circIndex(candArenaHeadings(mindex)-steps(k),nx);
                rCTemp = centerHeadingNeuronEB(theta,headingTemp);

                hUpdate  = hMat(arenaHeadingUpdate,:);
                dh       = relu((1-rCTemp)-hUpdate ).*nonlin(1-hUpdate ) - relu(hUpdate -(1-rCTemp)).*nonlin(hUpdate );
                hUpdate  = hUpdate  + alphaH.*dh.*abs( turnSize./turnDur).^2.*(turnDur./ddt);
                
                hMat(arenaHeadingUpdate,:) = hUpdate;
            end

        end
    end
    hMatAll = cat(3,hMatAll,hMat);

   
    
    %---------------------determine if bump will jump---------------------%
    if strcmp(envType,'sym')
        
        candArenaHeadings = [circIndex(arenaHeadingTemp-2:arenaHeadingTemp+2,nx),circIndex(arenaHeadingTemp+nx/2-2:arenaHeadingTemp+nx/2+2,nx)];
        activities = [.33,.66,1,.66,.33,.33,.66,1,.66,.33];
        
        dp1 = activities*hMat(candArenaHeadings,headingTemp)./numel(activities);
        dp2 = activities*hMat(candArenaHeadings,circIndex(headingTemp-nx/2,nx))./numel(activities);

        xx = [1./dp1,1./dp2];
        xx = xx./sum(xx);
        pJump(i) = 1-xx(1);

        if rand()<pJump(i)
            headingTemp = circIndex(headingTemp-nx/2,nx);
            jumps(i) = 1;
            jumpTraj(end) = 1;
            offsets = [offsets,offsets(end)];
        else
            jumps(i) = 0;
            offsets = [offsets,1-offsets(end)];
        end
    elseif strcmp(envType,'asym_train')
        offsets = [offsets,offsets(end)];
        
        %extract weights for current RN input
        candHeadings = circIndex(headingTemp-6:headingTemp+6,nx);
        xx = inf(1,nx);
        xx(candHeadings) = hMat(arenaHeadingTemp,candHeadings);
        xx = 1./xx;
        xx = xx./sum(xx);
        
        indheadings = 1:nx;
        %sort probabilities
        [ps,kk] = sort(xx);
        
        indselect = find(histcounts(rand(),[0,cumsum(ps)]));
        headingTemp = indheadings(kk(indselect));
        
        pJump(i) = nan;
        jumps(i) = nan;
    else
        offsets = [offsets,offsets(end)];
        pJump(i) = nan;
        jumps(i) = nan;
    end

    %--------------------------update heading-----------------------------%
    heading      = headingTemp;
    arenaHeading = arenaHeadingTemp;
    
  
end

res.weights.iter    = iterNum;
res.weights.goal    = wAllIter;
res.weights.heading = hMatAll;


res.behavior.arenaHeading = arenaHeadings;
res.behavior.heading      = headings;

res.behavior.bumpProfile = bumpProfileIter;
res.behavior.sacDir    = dirs;
res.behavior.sacProb   = pRight;
res.behavior.sacSize   = sizes;
res.behavior.fixDur    = durs;
res.behavior.bumpJump  = jumps;
res.behavior.pJump     = pJump;
res.behavior.pJumpBase = pJumpBase;
res.behavior.reward    = rewards;
res.behavior.safe      = safeVisits;

res.dtTot = dtTot;
res.headingBins = 1:(nx+1);
res.headingLabels = headingLabels;
res.arenaLabels   = arenaLabels;

res.behavior.headingTraj      = headingTraj;
res.behavior.arenaHeadingTraj = arenaHeadingTraj;
res.behavior.visSceneTraj     = visSceneTraj;
res.behavior.indsTraj         = indsTraj;
res.behavior.jumpTraj         = jumpTraj;

res.behavior.headingTrajU      = unwrap(headingTraj);
res.behavior.arenaHeadingTrajU = unwrap(arenaHeadingTraj);
res.behavior.visSceneTrajU     = unwrap(visSceneTraj);

res.behavior.offsets = offsets;


end

function xabs = unwrap(x)

%maps position onto absolute coordinates
xabs = x;
dx = diff(xabs);
ii=find(ismember(abs(dx),31));
for i=1:numel(ii)
    xabs( (ii(i)+1):end ) = xabs( (ii(i)+1):end ) - sign(dx(ii(i)))*32;
end

end

