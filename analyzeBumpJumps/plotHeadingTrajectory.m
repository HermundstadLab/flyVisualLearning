function plotHeadingTrajectory(tSeriesProcAll,fly)
% PLOTHEADINGTRAJECTORY plots the internal heading trajectory of the bump 
% in the EB using the fluorescence data from a single fly
%
% INPUTS:
%   tSeriesProcAll: data structure containing fluorescence data
%   fly:   integer index of fly to plot
%  

figure;set(gcf,'Position',[200 200 1400 600],'color','w')
subplot(3,1,1);
imagesc(tSeriesProcAll(fly).jumpData(1).dffMat);colormap(flipud(gray));set(gca,'ydir','normal');hold on;
plot(tSeriesProcAll(fly).jumpData(1).bumpjumps.locsMaxSinglePeaks./(2*pi/32)+16.5,'-k');set(gca,'ydir','normal');
set(gca,'fontsize',16)
title(['heatmap bounds (\Delta F/F): ', num2str(min(tSeriesProcAll(fly).jumpData(1).dffMat(:))),', ',num2str(max(tSeriesProcAll(fly).jumpData(1).dffMat(:)))])
ylabel('compass heading (wedges)')
xlabel('time (s)')
colorbar;


for trial = 1:2
    subplot(3,1,1+trial);
    
    x1 = tSeriesProcAll(fly).jumpData(trial).arena-tSeriesProcAll(fly).jumpData(trial).bumpjumps.offsets(1);
    x1(x1>pi) = x1(x1>pi)-2*pi;
    x2 = tSeriesProcAll(fly).jumpData(trial).arena-tSeriesProcAll(fly).jumpData(trial).bumpjumps.offsets(1)+pi;
    x2(x2>pi) = x2(x2>pi)-2*pi;
    plot(tSeriesProcAll(fly).jumpData(trial).bumpjumps.locsMaxSinglePeaks,'-k');hold on
    plot(x1)
    plot(x2)
    set(gca,'fontsize',16)
    ylabel('azimuth (deg)')
    xlabel('time (s)')
    xlim([0,numel(x1)])
    ylim([-pi,pi])
end
