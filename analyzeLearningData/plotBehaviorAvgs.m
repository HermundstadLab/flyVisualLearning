function [tFall,pSall,dxSall,vSall,fSall] = plotBehaviorAvgs(s,shift,align,p,pplot,windowsize,trials,additionalProperty)
% PLOTBEHAVIORAVGS plots the heading-dependent averages of different
% behavioral features, computed for a set of
%
% INPUTS:
%   s:  data structure containining segmented data
%   shift: string that specifies whether or not to collapse data to half 
%       of the arena. Can take values 'full', 'full_flat', or 'half'
%   align: string that specifies whether to align to the safe zone, or to
%       flies preferences. Can take values of 'pref' or 'safe'
%   p:  parameter vector
%   pplot: plotting parameter vector
%   windowsize: integer that specifies the width of bins to use for 
%       averaging. Default value is 7.
%   trials: array that specifies the trials to combine for averages
%   additionalProperty: string that specifies an additional behavioral
%       property to show alongside the duration of fixations and direction 
%       of saccades. Can take values of 'saccadeSize', 'saccadeVelocity',
%       'saccadeFrequency', or 'fixationFrequency'. Default is 
%       'saccadeSize'.    
%

if nargin<8
    additionalProperty = 'saccadeSize';
    if nargin<7
        trials = [p.naiveTrials, p.probeTrials1,p.probeTrials2];
        if nargin<6
            windowsize = 7;
        end   
    end
end

if strcmp(align,'pref')==1
    flies = 1:numel(s.fly);
elseif strcmp(align,'safe')==1
    for i=1:numel(s.fly)
        for j=1:numel(s.fly(1).trial)
            PI(i,j) = s.fly(i).trial(j).PI; 
        end
    end

    flies = 1:numel(s.fly);
else
    error('unrecognized alignment')
end
    
if strcmp(shift,'full')==1
    xmax = p.nx; 
elseif strcmp(shift,'full_flat')==1
    xmax = p.nx; 
elseif strcmp(shift,'half')==1
    xmax = p.nx/2;
else
    error('unrecognized shift')
end

tFall = [];
fFall = [];
pSall = [];
dxSallR = [];
dxSallL = [];
vSallR  = [];
vSallL  = [];
fSallR  = [];
fSallL  = [];
nFall = [];
nSall = [];

for i=1:numel(flies)
    [tFavg,fFavg,pSavg,dxSavgR,dxSavgL,vSavgR,vSavgL,fSavgR,fSavgL,nF,nS,bins] = getAvgs(s,i,shift,align,p,windowsize,trials);
    tFall = [tFall;tFavg];
    fFall = [fFall;fFavg];
    pSall = [pSall;pSavg];
    dxSallR = [dxSallR;dxSavgR];
    dxSallL = [dxSallL;dxSavgL];
    vSallR  = [vSallR;vSavgR];
    vSallL  = [vSallL;vSavgL];
    fSallR  = [fSallR;fSavgR];
    fSallL  = [fSallL;fSavgL];
    nFall = [nFall;nF];
    nSall = [nSall;nS];
end


%require a minimum number of samples to include
nsamplesmin = 2;
tFall(  nFall<nsamplesmin) = nan;
fFall(  nFall<nsamplesmin) = nan;
pSall(  nSall<nsamplesmin) = nan;
dxSallR(nSall<nsamplesmin) = nan;
dxSallL(nSall<nsamplesmin) = nan;
vSallR( nSall<nsamplesmin) = nan;
vSallL( nSall<nsamplesmin) = nan;
fSallR( nSall<nsamplesmin) = nan;
fSallL( nSall<nsamplesmin) = nan;


%require a minimum number of flies to include
nfliesmin = 5;
tFall(:,sum(~isnan(tFall),1)<nfliesmin) = nan;
fFall(:,sum(~isnan(fFall),1)<nfliesmin) = nan;
pSall(:,sum(~isnan(pSall),1)<nfliesmin) = nan;
dxSallR(:,sum(~isnan(dxSallR),1)<nfliesmin) = nan;
dxSallL(:,sum(~isnan(dxSallL),1)<nfliesmin) = nan;
vSallR( :,sum(~isnan(vSallR ),1)<nfliesmin) = nan;
vSallL( :,sum(~isnan(vSallL ),1)<nfliesmin) = nan;
fSallR( :,sum(~isnan(fSallL ),1)<nfliesmin) = nan;
fSallL( :,sum(~isnan(fSallL ),1)<nfliesmin) = nan;

dxSall = (dxSallR+dxSallL)./2;
vSall  = (vSallR +vSallL )./2; 
fSall  = (fSallR +fSallL )./2; 

cF = pplot.cF;
cS = pplot.cS;


    
figure;set(gcf,'Position',[200 200 800 1000],'color','w')


%------------------- plot duration of fixations ----------------------%
subplot(3,5,[1,2,3]);hold on;
y0 = tFall;

ya = nanmean(y0);
ys = nanstd(y0)./sqrt(sum(~isnan(y0),1));
plot(bins,y0,'linewidth',.5,'color',[.8,.8,.8])
plot(bins,ya,'linewidth',3,'color',cF)

ylim([0,30]);
yl = ylim;
xlim([0,xmax]);
plot([24,24],[yl(1),yl(2)],'--k')
plot([72,72],[yl(1),yl(2)],'--k')
xticks([0,24,48,72,96])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('heading')
ylabel('avg duration of fixations')
set(gca,'fontsize',16)

subplot(3,5,[4,5]);hold on;
if strcmp(shift,'full_flat')==1
    bins1 = bins(1:96);
    xx1 = [bins1,fliplr(bins1),bins1(1)];

    ya1 = ya(1:96);
    ys1 = ys(1:96);
    yy1 = [ya1-ys1,fliplr(ya1+ys1),ya1(1)-ys1(1)];

    fill(xx1,yy1,cF,'FaceAlpha',0.5,'edgecolor','none')
    plot(bins1,ya1,'linewidth',2,'color',cF)
    
    ylim([0,12])
    yl = ylim;
    yticks([0,4,8,12])
    plot([24,24],[yl(1),yl(2)],'--k')
    plot([72,72],[yl(1),yl(2)],'--k')
    xlim([0,96]);
    xticks([0,24,48,72,96])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
else
    
    bins1 = bins(1:48);
    xx1 = [bins1,fliplr(bins1),bins1(1)];
    ya1 = ya(1:48);
    ys1 = ys(1:48);
    yy1 = [ya1-ys1,fliplr(ya1+ys1),ya1(1)-ys1(1)];

    fill(xx1,yy1,cF,'FaceAlpha',0.5,'edgecolor','none')
    plot(bins1,ya1,'linewidth',2,'color',cF)

    ylim([0,12])
    yticks([0,4,8,12])
    xlim([0,48]);
    xticks([0,24,48])
    xticklabels({'0','\pi/2','\pi'})
    
    if strcmp(shift,'full')==1
        ya2 = ya(49:96);
        ys2 = ys(49:96);
        yy2 = [ya2-ys2,fliplr(ya2+ys2),ya2(1)-ys2(1)];

        fill(xx1,yy2,cF,'FaceAlpha',0.5,'edgecolor','none')
        plot(bins1,ya2,'linewidth',2,'color',cF)
    end
    yl = ylim;
    plot([24,24],[yl(1),yl(2)],'--k')
    xlim([0,48]);
    xticks([0,24,48])
    xticklabels({'0','\pi/2','\pi'})
end
xlabel('heading')
set(gca,'fontsize',16)

%------------------- plot bias of saccades ---------------------------%
subplot(3,5,[6,7,8]);hold on;
ya = nanmean(pSall);
ys = nanstd(pSall)./sqrt(sum(~isnan(pSall),1));
plot(bins,pSall,'linewidth',.5,'color',[.8,.8,.8])
plot(bins,ya,'linewidth',2,'color',cS)

ylim([-.75,.75]);
yl = ylim;
xlim([0,xmax]);
plot([24,24],[yl(1),yl(2)],'--k')
plot([72,72],[yl(1),yl(2)],'--k')
xticks([0,24,48,72,96])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('heading')
ylabel('probability of rightward saccade')
set(gca,'fontsize',16)


subplot(3,5,[9,10]);hold on;
if strcmp(shift,'full_flat')==1
    bins1 = bins(1:96);
    xx1 = [bins1,fliplr(bins1),bins1(1)];

    ya1 = ya(1:96);
    ys1 = ys(1:96);
    yy1 = [ya1-ys1,fliplr(ya1+ys1),ya1(1)-ys1(1)];

    fill(xx1,yy1,cS,'FaceAlpha',0.5,'edgecolor','none')
    plot(bins1,ya1,'linewidth',2,'color',cS)
    
    ylim([-.25,.25]);
    xlim([0,96]);
    xticks([0,24,48,72,96])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    yl = ylim;
    plot([24,24],[yl(1),yl(2)],'--k')
    plot([72,72],[yl(1),yl(2)],'--k')
    
else
    bins1 = bins(1:48);
    xx1 = [bins1,fliplr(bins1),bins1(1)];
    
    ya1 = ya(1:48);
    ys1 = ys(1:48);
    yy1 = [ya1-ys1,fliplr(ya1+ys1),ya1(1)-ys1(1)];

    fill(xx1,yy1,cS,'FaceAlpha',0.5,'edgecolor','none')
    plot(bins1,ya1,'linewidth',2,'color',cS)

    ylim([-.25,.25]);
    xlim([0,48]);
    xticks([0,24,48])
    xticklabels({'0','\pi/2','\pi'})
    yl = ylim;
    if strcmp(shift,'full')==1
        ya2 = ya(49:96);
        ys2 = ys(49:96);
        yy2 = [ya2-ys2,fliplr(ya2+ys2),ya2(1)-ys2(1)];

        fill(xx1,yy2,cS,'FaceAlpha',0.5,'edgecolor','none')
        plot(bins1,ya2,'linewidth',2,'color',cS)
        xticks([0,24,48])
    end
    plot([24,24],[yl(1),yl(2)],'--k')
end

xlabel('heading')
set(gca,'fontsize',16)



%--------- plot size, speed, freq of saccades or freq of fixations -------%
subplot(3,5,[11,12,13]);hold on;
if strcmp(additionalProperty,'saccadeSize')
    y0 = (dxSallR+dxSallL)./2;
    ylims = [.8,1.2];
elseif strcmp(additionalProperty,'saccadeVelocity')
    y0 = (vSallR+vSallL)./2;
    ylims = [.8,1.2];
elseif strcmp(additionalProperty,'saccadeFrequency')
    y0 = (fSallR+fSallL)./2;
    ylims = [.98,1.02];
elseif strcmp(additionalProperty,'fixationFrequency')
    y0 = fFall;
    ylims = [.98,1.02];
else
    error('unrecognized saccade property')
end

ya = nanmean(y0);
ys = nanstd(y0)./sqrt(sum(~isnan(y0),1));
plot(bins,y0,'linewidth',.5,'color',[.8,.8,.8])
plot(bins,ya,'linewidth',2,'color','k')

ylim([0.225,2.225]);
yl = ylim;
xlim([0,xmax]);
plot([24,24],[yl(1),yl(2)],'--k')
plot([72,72],[yl(1),yl(2)],'--k')
xticks([0,24,48,72,96])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('heading')
ylabel('net size (R+L)')
set(gca,'fontsize',16)

subplot(3,5,[14,15]);hold on;
if strcmp(shift,'full_flat')==1
    bins1 = bins(1:96);
    xx1 = [bins1,fliplr(bins1),bins1(1)];

    ya1 = ya(1:96);
    ys1 = ys(1:96);
    yy1 = [ya1-ys1,fliplr(ya1+ys1),ya1(1)-ys1(1)];

    fill(xx1,yy1,'k','FaceAlpha',0.5,'edgecolor','none')
    plot(bins1,ya1,'linewidth',2,'color','k')
    
    ylim(ylims);
    xticks([0,24,48,72,96])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    yl = ylim;
    plot([24,24],[yl(1),yl(2)],'--k')
    plot([72,72],[yl(1),yl(2)],'--k')
else
    bins1 = bins(1:48);
    xx1 = [bins1,fliplr(bins1),bins1(1)];

    ya1 = ya(1:48);
    ys1 = ys(1:48);
    yy1 = [ya1-ys1,fliplr(ya1+ys1),ya1(1)-ys1(1)];

    fill(xx1,yy1,'k','FaceAlpha',0.5,'edgecolor','none')
    plot(bins1,ya1,'linewidth',2,'color','k')

    ylim(ylims);
    xticks([0,24,48])
    xticklabels({'0','\pi/2','\pi'})
    
    if strcmp(shift,'full')==1
        ya2 = ya(49:96);
        ys2 = ys(49:96);
        yy2 = [ya2-ys2,fliplr(ya2+ys2),ya2(1)-ys2(1)];

        fill(xx1,yy2,'k','FaceAlpha',0.5,'edgecolor','none')
        plot(bins1,ya2,'linewidth',2,'color','k')
        ylim(ylims);
    end
    yl = ylim;
    plot([24,24],[yl(1),yl(2)],'--k')
end

xlabel('heading')
set(gca,'fontsize',16)

end





function [tFavg,fFavg,pSavg,dxSavgR,dxSavgL,vSavgR,vSavgL,fSavgR,fSavgL,nF,nS,bins] = getAvgs(s,fly,shift,align,p,ddx,trials)


if strcmp(shift,'half')==1
    nx = p.nx/2;
else
    nx = p.nx;
end
bins = 1:nx;

    
dtF = [];
fF  = [];
dxS = [];
vS  = [];
fS  = [];

%extract data
for j=1:numel(trials)
    trial = trials(j);
    if strcmp(align,'safe')==1
        if strcmp(shift,'half')==1
            sF = s.fly(fly).trial(trial).fixations.xavg;
            sS = s.fly(fly).trial(trial).saccades.xstart;
        else
            sF = s.fly(fly).trial(trial).fixations.xavgF;
            sS = s.fly(fly).trial(trial).saccades.xstartF;
        end
    else
        if strcmp(shift,'half')==1
            sF = s.fly(fly).trial(trial).shifted.fixations.xavg;
            sS = s.fly(fly).trial(trial).shifted.saccades.xstart;
        else
            sF = s.fly(fly).trial(trial).shiftedFull.fixations.xavgF;
            sS = s.fly(fly).trial(trial).shiftedFull.saccades.xstartF;
        end
    end

    if numel(sF)>0
        if all(~isnan(sF))
            if isfield(s.fly(fly).trial(trial).fixations,'dt')
                dt = s.fly(fly).trial(trial).fixations.dt./1000;
                df = s.fly(fly).trial(trial).fixations.wbfAvg;
                ii = find(dt>p.durThreshold);
                dtF  = [dtF;[sF(ii)',dt(ii)']];
                fF   = [fF; [sF(ii)',df(ii)']];
            end
            if isfield(s.fly(fly).trial(trial).saccades,'dt')
                dx = s.fly(fly).trial(trial).saccades.dx;
                dv = s.fly(fly).trial(trial).saccades.wbaAvg;
                df = s.fly(fly).trial(trial).saccades.dt./1000;
                jj = find(abs(dv)>p.velThreshold);
                dxS = [dxS;[sS(jj)',dx(jj)']];
                vS  = [vS; [sS(jj)',dv(jj)']];
                fS  = [fS; [sS(jj)',df(jj)']];
            end
        end
    end
end

%filter
if numel(dtF)>0
    for m = 1:nx
        inds = mod(m-ddx:m+ddx,nx);
        iF = find(ismember(mod(round(dtF(:,1)),nx),inds));
        tFavg(m) = mean(dtF(iF,2));
        fFavg(m) = mean(fF(iF,2))./mean(fF(:,2));
        nF(m) = numel(iF);
    end
else
    tFavg(1:nx) = nan;
    nF(1:nx) = nan;
end
if numel(dxS)>0
    for m=1:nx
        inds = mod(m-ddx:m+ddx,nx);
        iS  = find(ismember(dxS(:,1),inds));
        pSavg(m) = nanmean((sign(dxS(iS,2))+1)./2) - nanmean((sign(dxS(:,2))+1)./2);
        nS(m) = numel(iS);

        dxsR = dxS(iS(sign(dxS(iS,2))>0),2);
        dxsR = abs(dxsR)./nanmean(abs(dxS(sign(dxS(:,2))>0,2)));
        
        dxsL = dxS(iS(sign(dxS(iS,2))<0),2);
        dxsL = abs(dxsL)./nanmean(abs(dxS(sign(dxS(:,2))<0,2)));

        dxSavgR(m) = nanmean(dxsR);
        dxSavgL(m) = nanmean(dxsL);


        vsR = vS(iS(sign(vS(iS,2))<0),2);
        vsR = abs(vsR)./nanmean(abs(vS(sign(vS(:,2))<0,2)));
        
        vsL = vS(iS(sign(vS(iS,2))>0),2);
        vsL = abs(vsL)./nanmean(abs(vS(sign(vS(:,2))>0,2)));

        vSavgR(m) = nanmean(vsR);
        vSavgL(m) = nanmean(vsL);


        fsR = fS(iS(sign(vS(iS,2))<0),2);
        fsR = abs(fsR)./nanmean(abs(fS(sign(vS(:,2))<0,2)));
        
        fsL = fS(iS(sign(vS(iS,2))>0),2);
        fsL = abs(fsL)./nanmean(abs(fS(sign(vS(:,2))>0,2)));

        fSavgR(m) = nanmean(fsR);
        fSavgL(m) = nanmean(fsL);
       

    end
else
    pSavg(  1:nx) = nan;
    dxSavgR(1:nx) = nan;
    dxSavgL(1:nx) = nan;
    vSavgR( 1:nx) = nan;
    vSavgL( 1:nx) = nan;
    nS(1:nx) = nan;
end


end









