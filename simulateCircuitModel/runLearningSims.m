function [headLocs,goalLocs,headGains,goalGains,PIscores,distDanger,angdiff,headingVisits,headingJumps,res] = runLearningSims(p,gain,loc,params)
% RUNLEARNINGSIMS runs many simulations of an agent (model fly) that must form
% a mapping of the visual environment while learning to avoid danger zones.
% Simulations vary in the initial strength of the goal heading. 
%
% INPUTS:
%   p:      structure containing simulation parameters
%   gain:   initial strength of goal heading
%   loc:    initial location of goal heading
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


if nargin<4
    params = {'sym','train','rand_init'};
end
nx = p.N;
theta  = linspace(0,2*pi,nx+1);
theta  = theta(1:end-1);

headingLabels = 1:nx;
arenaLabels   = fliplr(headingLabels);

times = [1,5:5:1600]; 
nmax  = numel(times);

ICs   = [gain;loc];
niter = size(ICs,2);

nhist = p.N+1;

headLocs    = nan(niter,nmax);
goalLocs    = nan(niter,nmax);
headGains   = nan(niter,nmax);
goalGains   = nan(niter,nmax);
PIscores    = nan(niter,nmax);   
distDanger  = nan(niter,nmax); 
angdiff     = nan(niter,nmax); 
    
headingJumps  = nan(niter,nhist);
headingVisits = nan(niter,nhist);

for m=1:niter

    %initialize heading weights in asymmetric environment

    %strength of goal vector (min val = 0, max val = .5)
    goalStrength = ICs(1,1);

    %fix goal weights
    indGoal = ICs(2,1);

    [~,indGoal_heading] = find(headingLabels==indGoal);
    r.w = goalStrength*cos(theta-theta(indGoal_heading))+.5;


    gMat  = r.w;
    rW    = sum(gMat.*exp(1i.*theta));
    angW  = angle(rW);
    goalW = mod(angW,2*pi); 
    ss = abs(rW)./12;
    dd = min(abs(goalW - [pi/2,3*pi/2]));

    %randomly initialize heading weights
    r.hMat = .25*rand(nx)+.5;

    %initial headings
    r.heading = randperm(nx,1);  

    %trial type
    r.trialType = 'probe';
    r.envType   = 'asym';

    r.dtTot = 1000;
    r.mmax  = 10000;

    %run single simulation to initialize weights
    if strcmp(params{3},'asym_init')
        r.envType   = 'asym';
        res    = simulateCircuit(p,r);
        r.hMat = res.weights.heading(:,:,end);
    elseif strcmp(params{3},'sym_init')
        r.envType   = 'sym';
        res    = simulateCircuit(p,r);
        r.hMat = res.weights.heading(:,:,end);
    end


    %run second simulation in symmetric environment
    %initial heading
    r.heading = randperm(nx,1);  

    %trial type
    r.envType   = params{1};
    r.trialType = params{2};
    
    r.dtTot = 2000;
    r.mmax  = 20000;

    res = simulateCircuit(p,r);

    %compute empirical histogram of bump jumps (as written, only useful for fixed goal heading)
    heading = res.behavior.heading;
    headingJumps(m,:)  = histcounts(heading(res.behavior.bumpJump>.5),0:nhist);
    headingVisits(m,:) = histcounts(heading,0:nhist);


    for j=1:numel(times)+1


        if j==1
            ind = 1;
        else
            ind = find(cumsum(res.behavior.fixDur)>times(j-1),1,'first');
        end

        gMat  = res.weights.goal(ind,:);
        rW    = sum(gMat.*exp(1i.*theta));
        angW  = angle(rW);
        goalW = mod(angW,2*pi); 
        strengthW = abs(rW)./12;

        hMat = squeeze(res.weights.heading(:,:,ind));

        pJ = [];
        for k = 1:nx
            %k = internal heading
            %h = arena heading:
            h = arenaLabels(k);
            candArenaHeadings = [circIndex(h-1:h+1,nx),circIndex(h+nx/2-1:h+nx/2+1,nx)];
            activities = [.5,1,.5,.5,1,.5];

            dp1 = activities*hMat(candArenaHeadings,k)./numel(activities);
            dp2 = activities*hMat(candArenaHeadings,circIndex(k-nx/2,nx))./numel(activities);

            %compute pStay, not pJump
            xx = [1./dp1,1./dp2];
            xx = xx./sum(xx);
            pJ(k) = xx(1); %this is pStay

        end

        rH    = sum(pJ.*exp(1i.*theta));
        angH  = angle(rH);
        goalH = mod(angH,2*pi);
        strengthH = abs(rH)./12;

        goalLocs(m,j)  = goalW;
        headLocs(m,j)  = goalH;

        distDanger(m,j) = min(abs(goalW - [pi/2,3*pi/2]));

        goalGains(m,j) = strengthW;
        headGains(m,j) = strengthH;

        hvec = [strengthH.*cos(angH), strengthH.*sin(angH)];
        rvec = [strengthW.*cos(angW), strengthW.*sin(angW)];
        angdiff(m,j)   = acos(sum(hvec.*rvec)./(strengthW.*strengthH));


        %compute PI score
        %initialize goal weights
        rtmp.w = gMat;

        %initialize heading weights
        rtmp.hMat = hMat;

        %initial headings
        rtmp.heading = randperm(nx,1);  

        %trial type
        rtmp.trialType = 'frozen';
        rtmp.envType   = params{1};

        rtmp.dtTot = 1000;
        rtmp.mmax  = 10000;

        rout = simulateCircuit(p,rtmp);
        num = sum(rout.behavior.fixDur(rout.behavior.safe>0)) - sum(rout.behavior.fixDur(rout.behavior.safe<0));
        den = sum(rout.behavior.fixDur);
        PIscores(m,j) = num./den;

    end
end


end


