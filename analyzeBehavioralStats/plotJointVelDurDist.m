function plotJointVelDurDist(sL,sNL,p)
% PLOTJOINTVELDURDIST plots the joint distribution of angular velocities
% and durations for saccades and fixations
%
% INPUTS:
%   sL, sNL:  data structures containining segmented data
%   p: structure containing parameter values
% 

datasets = {sL,sNL};
nSets    = numel(datasets);
nFlies   = [numel(sL.fly),numel(sNL.fly)];
nTrials  = [p.nTrials,p.nTrials];

durThreshold = p.durThreshold;
velThreshold = p.velThreshold;

fixV  = getDists(datasets,'trials','fixations','wbaAvg',nTrials,nFlies);
fixT  = getDists(datasets,'trials','fixations','dt',nTrials,nFlies);
sacV  = getDists(datasets,'trials','saccades','wbaAvg',nTrials,nFlies);
sacT  = getDists(datasets,'trials','saccades','dt',nTrials,nFlies);
sacX  = getDists(datasets,'trials','saccades','dx',nTrials,nFlies);

FV = [];
FT = [];
SV = [];
ST = [];
SX = [];
for i=1:nSets
    for j=1:nTrials(i)
        FV = [FV;fixV{i,j}'];
        FT = [FT;fixT{i,j}'];
        SV = [SV;sacV{i,j}'];
        ST = [ST;sacT{i,j}'];
        SX = [SX;sacX{i,j}'];
    end
end
ii = find(FT<durThreshold);
jj = find(abs(SV)<velThreshold);

FV(ii) = [];
FT(ii) = [];
SV(jj) = [];
ST(jj) = [];
SX(jj) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              plot joint distribution of velocity and duration           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
cF = [0,160,126]./255;
cS = [175,175,255]./255;
subplot(2,3,1);hold on;

%coordinates for zoomed-in region
x = [-2, 2, 2, -2, -2];
y = [ 0, 0, 5,  5,  0];
scatter(FV,FT,100,cF,'filled')
scatter(SV,ST,100,cS,'filled')
plot(x,y,'--k')
ylabel('duration \Delta T (s)')
xlabel('angular velocity \Delta WBA')
set(gca,'fontsize',16);


subplot(2,3,2);hold on;
scatter(FV,FT,100,cF,'filled')
xlim([-2,2])
ylim([0,5])
ylabel('duration \Delta T (s)')
xlabel('angular velocity \Delta WBA')
set(gca,'fontsize',16);

subplot(2,3,3);hold on;
scatter(SV,ST,100,cS,'filled')
xlim([-2,2])
ylim([0,5])
ylabel('duration \Delta T (s)')
xlabel('angular velocity \Delta WBA')
set(gca,'fontsize',16);

load('fits/sacDur.mat','vel')
nv = 16;
v = FV;
n = numel(v);
pvec = (0:(n-1))./(n-1);
vs = sort(v);
[~,i1] = min(abs(pvec-.025));
[~,i2] = min(abs(pvec-.975));

subplot(2,3,4);hold on;
f = (1:nv);
f = f./max(f);
for i=1:nv
    fill([vel(i),vel(i+1),vel(i+1),vel(i)],[0,0,1,1],f(i).*cS)
    fill([-vel(i),-vel(i+1),-vel(i+1),-vel(i)],[0,0,1,1],f(i).*cS)
end
for i=nv+1:numel(vel)-1
    fill([vel(i),vel(i+1),vel(i+1),vel(i)],[0,0,1,1],f(end).*cS)
    fill([-vel(i),-vel(i+1),-vel(i+1),-vel(i)],[0,0,1,1],f(end).*cS)
end

fill([vs(i1),vs(i2),vs(i2),vs(i1)],[1.5,1.5,2.5,2.5],cF)
xlim([-2, 2])
ylim([-.5,3])
set(gca,'fontsize',16);


d1 = fitdist(abs(SX),'lognormal');
d2 = fitdist(abs(SX*360/p.nx),'lognormal');

dx1med = exp(d1.mu);
dx2med = exp(d2.mu);

subplot(2,3,5);hold on;
scatter(-50.*SV.*ST,SX,100,cS,'filled')
plot([-250,250],[-250,250],'--k')
xlim([-250,250])
ylim([-250,250])
ylabel('\Delta x (pixels)')
xlabel('-50*\Delta T (s)*\Delta WBA')
title({'lognormal fits to \Delta x distribution:';...
    ['pixels: \mu = ',num2str(round(d1.mu,2)),', \sigma = ',num2str(round(d1.sigma,2)),...
    ' ; median size = ',num2str(round(dx1med))];...
    ['degrees: \mu = ',num2str(round(d2.mu,2)),', \sigma = ',num2str(round(d2.sigma,2)),...
    ' ; median size = ',num2str(round(dx2med,2))]})
set(gca,'fontsize',16);



set(gcf,'renderer','Painters')
set(gcf,'Position',[200 200 1600 800],'color','w')
    
    

    