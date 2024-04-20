function plotOptimalPolicy(p,pplot,policyType)
% PLOTOPTIMALPOLICY plots the final policy of an optimal agent trained to
% maintain a single goal heading.
%
% INPUTS:
%   p: structure containing simulation parameters
%   pplot: structure containing plotting parameters
%   policyType: string that specifies whether the policy should be trained
%       on an internal representation of heading, or directly on the visual
%       scene. Can take value 'heading' or 'visual'; default is 'heading'
%

if nargin<2
    policyType = 'heading';
end

p.nx = p.nx/2;

theta = linspace(0,2*pi,p.nx+1);
theta = theta(1:end-1);
offset = pi;
gainB  = p.gainB;
shiftB = p.shiftB;

durF = [];
dirS = [];
driftF  = [];
durFeff = [];
dirSeff = [];
nIter = 100;
parfor i=1:nIter
    [omega,fk,params] = trainFullPolicy(p,policyType);
    [nuF,pR,~,~] = convertParams(omega,fk,params);
    driftF  = [driftF;nuF];
    durF    = [durF;p.aF./nuF];
    durFeff = [durFeff;p.aF./(nuF.*(1-pBumpJump(theta,offset,gainB,shiftB)) + circshift(nuF,p.nx/2).*pBumpJump(theta,offset,gainB,shiftB))];
    
    dirS = [dirS;2.*pR-1];
    dirSeff = [dirSeff;2.*(pR.*(1-pBumpJump(theta,offset,gainB,shiftB)) + circshift(pR,p.nx/2).*pBumpJump(theta,offset,gainB,shiftB)) - 1];
end



%average direction and duration
driftFavg = mean(driftF);
durFavg   = mean(durF);
dirSavg   = mean(dirS);

durFjump = mean(durFeff);
dirSjump = mean(dirSeff);


if strcmp(policyType,'heading')==1
    %shift so that preferred heading is at different orientations (for visual comparison to data)
    theta = linspace(0,2*pi,p.nx+1);
    shiftComp   = 0;
    shiftArena  = p.nx/4;
    driftFavg0  = [driftFavg(shiftComp+1:end),driftFavg(1:shiftComp+1)];
    durFavg0    = [durFavg(shiftComp+1:end),durFavg(1:shiftComp+1)];
    dirSavg0    = [dirSavg(shiftComp+1:end),dirSavg(1:shiftComp+1)];
    durFavgJump = [durFjump(shiftArena+1:end),durFjump(1:shiftArena+1)];
    dirSavgJump = [dirSjump(shiftArena+1:end),dirSjump(1:shiftArena+1)];
    
    prefCompHead  = (p.nx/2-shiftComp )*2*pi/p.nx;
    prefArenaHead = (p.nx/2-shiftArena)*2*pi/p.nx;
    figure;set(gcf,'Position',[200 200 1400 800],'color','w');hold on;
    subplot(3,3,1);hold on;
    bar(theta,driftFavg0-mean(driftFavg0),'edgecolor','none','facecolor',pplot.cF)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('drift rate')
    xlabel('compass heading')
    yl = ylim;
    plot([prefCompHead,prefCompHead],yl,'--k')
    set(gca,'fontsize',16);
    
    subplot(3,3,2);hold on;
    plot(theta,durFavg0,'linewidth',pplot.lw,'color',pplot.cF)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('average duration of fixations')
    xlabel('compass heading')
    yl = ylim;
    plot([prefCompHead,prefCompHead],yl,'--k')
    set(gca,'fontsize',16);
    
    subplot(3,3,8);hold on;
    plot(theta(1:p.nx/2),durFavg0(1:p.nx/2),'linewidth',pplot.lw,'color',.5*pplot.cF)
    plot(theta(1:p.nx/2),durFavg0(p.nx/2+1:p.nx),'linewidth',pplot.lw,'color',pplot.cF)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('average duration of fixations')
    xlabel('arena heading')
    yl = ylim;
    plot([prefArenaHead,prefArenaHead],yl,'--k')
    set(gca,'fontsize',16);
    

    subplot(3,3,3);hold on;
    plot(theta,durFavgJump,'linewidth',pplot.lw,'color',pplot.cF)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('average duration of fixations')
    xlabel('arena heading')
    yl = ylim;
    plot([prefArenaHead,prefArenaHead],yl,'--k')
    set(gca,'fontsize',16);

    
    subplot(3,3,4);hold on;
    bar(theta,-dirSavg0+mean(dirSavg0),'edgecolor','none','facecolor',pplot.cS)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    yl = ylim;
    plot([prefCompHead,prefCompHead],yl,'--k')
    xlabel('compass heading')
    ylabel('directional bias of saccades')
    set(gca,'fontsize',16);
    
    subplot(3,3,5);hold on;
    plot(theta,dirSavg0,'linewidth',pplot.lw,'color',pplot.cS)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('directional bias of saccades')
    xlabel('compass heading')
    yl = ylim;
    plot([prefCompHead,prefCompHead],yl,'--k')
    set(gca,'fontsize',16);
    
    subplot(3,3,9);hold on;
    plot(theta(1:p.nx/2),dirSavg0(1:p.nx/2),'linewidth',pplot.lw,'color',.5*pplot.cS)
    plot(theta(1:p.nx/2),dirSavg0(p.nx/2+1:p.nx),'linewidth',pplot.lw,'color',pplot.cS)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('directional bias of saccades')
    xlabel('arena heading')
    yl = ylim;
    plot([prefArenaHead,prefArenaHead],yl,'--k')
    set(gca,'fontsize',16);
    
    
    subplot(3,3,6);hold on;
    plot(theta,dirSavgJump,'linewidth',pplot.lw,'color',pplot.cS)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('directional bias of saccades')
    xlabel('arena heading')
    yl = ylim;
    plot([prefArenaHead,prefArenaHead],yl,'--k')
    set(gca,'fontsize',16);

    
    subplot(3,3,7);hold on;
    plot(theta,pBumpJump(theta,offset+shiftComp,gainB,shiftB),'linewidth',pplot.lw,'color',pplot.cB)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    yl = ylim;
    plot([prefCompHead,prefCompHead],yl,'--k')
    xlabel('compass heading')
    ylabel('probability of bump jump')
    set(gca,'fontsize',16);
   
    
    
    

elseif strcmp(policyType,'visual')==1
    theta = linspace(0,2*pi,2*p.nx+1);
    durFavg = [durFavg,durFavg,durFavg(1)];
    dirSavg = [dirSavg,dirSavg,dirSavg(1)];
    
    figure('color','white','Units', 'Normalized', 'OuterPosition', [.2, 0.2, 0.2, .55]); hold on;
    subplot(2,1,1);hold on;
    plot(theta,durFavg,'linewidth',pplot.lw,'color',pplot.cF)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    ylabel('average duration of fixations')
    yl = ylim;
    plot([pi/2,pi/2],yl,'--k')
    set(gca,'fontsize',16);

    subplot(2,1,2);hold on;
    plot(theta,dirSavg,'linewidth',pplot.lw,'color',pplot.cS)
    xticks([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    yl = ylim;
    plot([pi/2,pi/2],yl,'--k')
    xlabel('heading')
    ylabel('directional bias of saccades')
    set(gca,'fontsize',16);

else
    error('unrecognized policy type')
end



