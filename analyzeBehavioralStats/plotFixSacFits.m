function plotFixSacFits(sL,p)
% PLOTFIXSACFITS plots the best-fitting parametric distributions to the
% distributions of saccade and fixation properties
%
% INPUTS:
%   sL:  data structures containining segmented data
%   p: structure containing parameter values
%

datasets = {sL};
nFlies   = numel(sL.fly);
nTrials  = p.nTrials;
durThreshold = p.durThreshold;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 plot fit params for fixation durations                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
figure;
subplot(2,2,1);hold on;

load('fits/fixDurTrials.mat')
ma = {'o','s'};

%extract reflexive fixations
fixT   = getDists(datasets,'trials','fixations','dt',nTrials,nFlies);
fixTD  = getDists(datasets,'trials','fixations','dtDanger',nTrials,nFlies);
fixTS  = getDists(datasets,'trials','fixations','tstart',nTrials,nFlies);
     
T  = fixT{ 1,3}';
TD = fixTD{1,3}';
TS = fixTS{1,3}';

ii  = find(TS<60000 & TD./T>.95 & T>durThreshold);
pdL = fitdist(T(ii),'inverse Gaussian');


% plot empirical means and variances
x = log(empMean.^3);
y = log(empVar);
for i=1:2
    scatter(x(i,:),y(i,:),200,'filled',ma{i},...
        'markerfacecolor','w','markeredgecolor',[.8,.8,.8],'linewidth',.75)
end

% plot fit parameters 
x = log(muT.^3);
y = log(muT.^3./lambdaT);
for i=1:2
    scatter(x(i,:),y(i,:),200,'filled',ma{i},...
        'markerfacecolor','w','markeredgecolor',[0,160,126]./255,'linewidth',2)
end

xL = log(mean(T(ii)).^3);
yL = log(var(T(ii)));
scatter(xL,yL,300,'filled','p','markerfacecolor','w','markeredgecolor',[.8,.8,.8],'linewidth',.75)

xL = log(pdL.mu.^3);
yL = log(pdL.mu.^3./pdL.lambda);
scatter(xL,yL,300,'filled','p','markerfacecolor','w','markeredgecolor','r','linewidth',2)



%find the equation of a line for variable nu, constant a
x0 = 2:.01:6;
b0 = mean(y(:)-x(:));
y0 = x0+b0;

%find the equation of a line for variable a, constant nu
a0 = exp(-b0/2);
muT0 = 4.5;
nu0 = a0./(exp(muT0).^(1/3));

%a = (a0-1):.01:(a0+1);
a = max(0.01,(a0-2)):.01:(a0+2);
s = 1;
x1 = log((a./nu0).^3);
y1 = log((a./nu0).^3.*s.^2./(a.^2));

plot(x0,y0,'--k')
plot(x1,y1,'--k')
xlim([3.5,6])
ylim([2,5])

lambdaT = -b0;
eta = 1;
a  = sqrt(eta.^2.*lambdaT);

nuR = a./pdL.mu;

xlabel('log(mean^3)')
ylabel('log(variance)')
set(gca,'fontsize',16)
title({'fit params for fixation duration';...
    ['\lambda = ',num2str(round(lambdaT,2)),...
    '; a = ',num2str(round(a,2)), ', \eta = ',num2str(eta),', \nu_R = ',num2str(round(nuR,2))]})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 plot fit params for saccade durations                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------ plot fit params for duration of saccades ---------------%
    

subplot(2,2,2);hold on
load('fits/sacDur.mat')

nv = 16;

% plot empirical means / variances
x = log(empMean.^3);
y = log(empVar);
for i=1:nv
    scatter(x(i),y(i),200,'filled',...
        'markerfacecolor','w','markeredgecolor',[.8,.8,.8],'linewidth',.75)
end

% plot fit params
x = log(muT.^3);
y = log(muT.^3./lambdaT);
f = (1:nv);
f = f./max(f);
for i=1:nv
    scatter(x(i),y(i),200,'filled',...
        'markerfacecolor',f(i).*[175,175,255]./255,'markeredgecolor',f(i).*[175,175,255]./255,'linewidth',2)
end

%find the equation of a line for variable nu, constant a, sig
%x0 = -3.6:.01:1.6;
x0 = -4:.01:1;
b0 = mean(y-x);
y0 = x0+b0;

%find the equation of a line for variable a, constant nu
a0 = exp(-b0/2);
muT0 = -3.15;
nu0 = a0./(exp(muT0).^(1/3));

%a = (a0-1):.01:(a0+1);
a = max(0.01,(a0-2)):.01:(a0+2);
s = 1;
x1 = log((a./nu0).^3);
y1 = log((a./nu0).^3.*s.^2./(a.^2));

lambdaT = -b0;
eta = 1;
a  = sqrt(eta.^2.*lambdaT);

plot(x0,y0,'--k')
plot(x1,y1,'--k')

ylim([-4.5,-1.4])
xlim([-4,-1])
xlabel('log(mean^3)')
ylabel('log(variance)')
set(gca,'fontsize',16);
title({'fit params for saccade duration';...
    ['\lambda = ',num2str(round(lambdaT,2)),...
    '; a = ',num2str(round(a,2)), ', \eta = ',num2str(eta)]})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                plot drift velocity vs saccade velocity                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



subplot(2,2,4);hold on;
load('fits/sacDur.mat');
fun = sigmoidFit;
v0 = 0:.01:2;
for i=1:nv
    scatter(vel(i),nu(i),200,'filled',...
        'markerfacecolor',f(i).*[175,175,255]./255,'markeredgecolor',f(i).*[175,175,255]./255,'linewidth',2)
end
for i=nv+1:numel(vel)
    scatter(vel(i),nu(i),200,'filled',...
        'markerfacecolor',f(end).*[175,175,255]./255,'markeredgecolor',f(end).*[175,175,255]./255,'linewidth',2)
end
plot(v0,fun(xlsq,v0),'--k','linewidth',1.5)

xlabel('saccade velocity v')
ylabel('drift velocity \nu' )
set(gca,'fontsize',16);
title({'drift velocity vs saccade velocity';...
    ['f_M = ',num2str(round(xlsq(1),2)),...
    ', k = ',num2str(round(xlsq(2),2)),...
    ', \omega_0 = ',num2str(round(xlsq(3),2)),...
    ', f_0 = ',num2str(round(xlsq(4),2))]})    
    
set(gcf,'Position',[200 200 1200 1000],'color','w')