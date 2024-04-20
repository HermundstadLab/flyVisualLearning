function plotDoublyPeakedTuningCurves(tSeriesProcAll,fly,trials)
%PLOTDOUBLYPEAKEDTUNINGCURVES plots tuning curves for EPG neurons in a
%symmetric visual scene
%
% INPUTS:
%   tSeriesProcAll: data structure containing fluorescence data
%   fly:   integer index of fly to plot
%   trials:  integer indices of trials to use
%  

M = [];
for i=trials
    M = cat(3,M,tSeriesProcAll(fly).jumpData(i).FRMaps.FRMapHeading);
end

FRmap = fliplr(rot90(squeeze(mean(M,3))));

figure;set(gcf,'Position',[200 200 800 1200],'color','w')
subplot(4,1,[1,2,3]);
imagesc(FRmap);colormap(flipud(gray));set(gca,'fontsize',16')
title(['heatmap bounds (\Delta F/F): ', num2str(min(FRmap(:))),', ',num2str(max(FRmap(:)))])
colorbar;
xlabel('arena heading (pixels)')
ylabel('compass heading (wedges)')

inds = 8:4:24;
cmap = gray(numel(inds)+2);
cmap = [cmap(2:3,:);cmap(1,:);cmap(4:5,:)]; 
for i=1:numel(inds)
    subplot(4,1,4);hold on;set(gca,'fontsize',16');xlim([0,48])
    if i==(numel(inds)+1)/2
        plot(FRmap(inds(i),:),'-','color',cmap(i,:),'linewidth',4)
    else
        plot(FRmap(inds(i),:),'-','color',cmap(i,:),'linewidth',2)
    end
    xlabel('arena heading (pixels)')
    ylabel('\Delta F/F')

    subplot(4,1,[1,2,3]);hold on;
    plot(1,inds(i),'x','color',cmap(i,:),'linewidth',20)
end
set(gcf,'color','w')
