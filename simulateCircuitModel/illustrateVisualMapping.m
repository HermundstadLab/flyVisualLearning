function illustrateVisualMapping(p,pplot,plotNum,params)
% ILLUSTRATEVISUALMAPPING generates different figures that schematize the
% visual mapping from ring neurons onto compass neurons:
%
% INPUTS
%   p: structure containing simulation parameters
%   pplot:  structure containing plotting parameters
%   plotNum key:
%       1:  external readout of residency, fixations, and saccades after
%           training heading weights
%       2:  simulation of trajectory with jumps
%   params: cell array of parameters that specify initial conditions
%       % params{1}: 'asym' or 'sym' (whether current scene is asymmetric or symmetric)
%       % params{2}: integer index that specifies location of goal heading
%       % params{3}: numeric value in [0,0.5] that specifies strength of goal heading
%       

if nargin<4
    params = {'asym',24,.5};
end

nx = p.N;
thetaF  = linspace(0,2*pi,nx+1);
theta  = thetaF(1:end-1);

headingLabels = 1:nx;
arenaLabels   = fliplr(headingLabels);
[~,isort] = sort(arenaLabels); 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,plotNum)  
    
    %params{1}: environment type ('sym', 'symNoJump', 'asym', of 'asym_train')
    %params{2}: goal location (either 8 or 24; center of safe zone)
    %params{3}: strength of goal vector (min val = 0, max val = .5)

    goalStrength = params{3};
    
    %fix goal weights
    indGoal =  params{2};
    [~,indGoal_heading] = find(headingLabels==indGoal);
    [~,indGoal_arena  ] = find(arenaLabels  ==indGoal);
    r.w = (goalStrength*cos(theta-theta(indGoal_heading))+.5);

    %randomly initialize heading weights
    r.hMat = .25*rand(nx)+.5;

    %initial headings
    r.heading = randperm(nx,1);  

    %trial type
    r.trialType = 'probe';
    r.envType   = params{1};

    r.dtTot = 5000;
    r.mmax  = 50000;

    
    cmax = 1;
    cmin = 0;
    
    %run single simulation to initialize weights
    res = simulateCircuit(p,r);
    
    nT = size(res.weights.heading,3);
    figure;set(gcf,'Position',[200 200 600 500],'color','w')

    %plot final weights
    M = squeeze(res.weights.heading(:,:,end));
    M = circshift(M,[0,1]);
    M = [M,M(:,1)];
    M = [M;M(1,:)];
    imagesc(theta,theta,M);caxis([cmin,cmax]);colormap(gray);hold on;
    plot([theta(indGoal_heading),theta(indGoal_heading)],[0,2*pi],'--k','linewidth',2)

    xlabel('compass heading')
    ylabel('arena heading')
    set(gca,'fontsize',16,'ydir','normal');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,plotNum)  
    
    cR = [174,0,66]./255;
    cP = [224,159,185]./255;


    %randomly initialize heading weights
    r.hMat = .25*rand(nx)+.5;

    %initial headings
    r.heading = randperm(nx,1);  

    %trial type
    r.envType   = 'sym';
    r.trialType = 'probe';

    r.dtTot = 5000;
    r.mmax  = 50000;
    

    indGoal = 8+nx/2; %center of safe zone;
    [~,indGoal_heading] = find(headingLabels==indGoal);
    [~,indGoal_arena  ] = find(arenaLabels  ==indGoal);
    
    r.w = (.5*cos(theta-theta(indGoal_heading))+.5);
    r.hMat    = .25*rand(nx)+.5;
    r.heading = randperm(nx,1);
    
    
    res = simulateCircuit(p,r);
    M = squeeze(res.weights.heading(:,:,end));
    r.hMat = M;
    
    %run another simulation to plot bump jumps
    res = simulateCircuit(p,r);
    
    figure;set(gcf,'Position',[200 200 1600 300],'color','w')
    tmax = 2000;
    t = 1:tmax;

    subplot(1,3,[1,2]);hold on;
    kk = find(res.behavior.jumpTraj(t)>0);
    for i=1:numel(kk)
        plot(t(kk(i):kk(i)+1),res.behavior.headingTraj(kk(i):kk(i)+1),'-k')
    end
    xlabel('time')
    ylabel('compass heading (wedges)')

    visScene180 = mod((res.behavior.visSceneTraj(t)-1)+16,32)+1;
    plot(t,res.behavior.visSceneTraj(t),'o','markerfacecolor',cR,'markeredgecolor','none');
    plot(t,visScene180,'o','markerfacecolor',cP,'markeredgecolor','none');
    plot(t,res.behavior.headingTraj(t), 'ok','markerfacecolor','k', 'markeredgecolor','none');
    set(gca,'fontsize',16)
    ylim([0,32])

    
    %grab final weights to use as initial values
    %for new sims (to read out the average behavior given this set of
    %weights, averaged over many short simulations)
    r.dtTot = 40;
    r.mmax  = 400;
    r.trialType = 'frozen';


    for i=1:1000

        r.heading = randperm(nx,1);  
        res = simulateCircuit(p,r);


        dd(1) = mean(res.behavior.offsets);
        dd(2) = 1-dd(1);
        d(i,1:2) = sort(dd,'descend');

    end

    subplot(1,3,3);hold on;
    d = mean(d);
    
    offsets(1) = thetaF(indGoal_arena);
    offsets(2) = thetaF(circIndex(indGoal_arena+nx/2,nx));
    r = 1;
    xx = linspace(0,2*pi,100);
    plot(r*cos(xx),r*sin(xx),'linewidth',.5,'color',[.8,.8,.8])
    plot(r/2*cos(xx),r/2*sin(xx),'linewidth',.5,'color',[.8,.8,.8])
    plot(r.*[-1,1],[0,0],'linewidth',.5,'color',[.8,.8,.8])
    plot([0,0],r.*[-1,1],'linewidth',.5,'color',[.8,.8,.8])
    for j=1:numel(d)
        plot([0,d(j).*cos(offsets(j))],[0,d(j).*sin(offsets(j))],'-k','linewidth',2)

        %set color = 1-d to match above colormap
        scatter(d(j).*cos(offsets(j)),d(j).*sin(offsets(j)),400*d(j),1-d(j),'filled');
        caxis([0,1]);colormap(gray)
    end
    pbaspect([1 1 1])
    axis off
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
