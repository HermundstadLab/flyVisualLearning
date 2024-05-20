function plotSims_VisualMapping(p,headGains,headingJumps,headingVisits,PIscores,gain,plotMaps,params)
% PLOTSIMS_VISUALMAPPING plots the output of simulations designed to probe 
% the evolution of the visual mapping from the visual scene onto the 
% compass heading
%
% INPUTS:
%   headGains:      strength of most stable compass heading over time
%   headingVisits:  histogram of visits to different compass headings
%   headingJumps:   histogram of headings at which compass heading jumped
%   PIscores:       PI score over time
%   gain:           array containing initial strengths of goal heading
%   plotMaps:       logical specifying whetehr to plot the heading map.
%                       Default is false.
%   params:         cell array specifying initial conditions
%       params{1}: 'asym' or 'sym' (whether current scene is asymmetric or symmetric)
%       params{2}: 'train' or 'probe' (whether or not model flies experience punishment)
%       params{3}: 'sym_init' or 'asym_init' (whether visual map is initialized in a symmetric or asymmetric scene)
%

if nargin<7
    plotMaps = false;
    params = {'sym','probe','asym_init'};
end

gain = reshape(gain,[200*200,5]);
gain = gain(1:200,:);
gain = gain(:);
gainU = unique(gain);

cmap = flipud(gray(6));
cmap = cmap(2:end,:);

figure;hold on;set(gcf,'Position',[200 200 1200 600],'color','w')
subplot(2,3,1);hold on;
for i=1:5
    ii = find(gain==gainU(i));
    num = headingJumps(ii,2:end);den = headingVisits(ii,2:end);
    num = [num(:,1),num(:,2:16)+fliplr(num(:,18:32)),num(:,17)];
    den = [den(:,1),den(:,2:16)+fliplr(den(:,18:32)),den(:,17)];
    plot(0:16,fliplr(nanmean(num./den)),'linewidth',2,'color',cmap(i,:));
end
xlim([0,16])
xticks([0,8,16])
xticklabels({'0','90','180'})
set(gca,'fontsize',16)
xlabel('distance to goal heading (deg)')
ylabel('conditional probability of jump')



for i=1:5

    subplot(2,3,2);hold on;
    set(gca,'fontsize',16)
    xlabel('final strength of visual map')
    ylabel('final strength of behavioral preference')
    ii = find(gain==gainU(i));
    h = headGains(ii,:);h = h(:);
    pp = -PIscores(ii,:);pp = pp(:);
    plot(mean(h),mean(pp),'o','markeredgecolor','none','markerfacecolor',cmap(i,:));
    plot([mean(h)-std(h),mean(h)+std(h)],[mean(pp),mean(pp)],'-','color',cmap(i,:))
    plot([mean(h),mean(h)],[mean(pp)-std(pp),mean(pp)+std(pp)],'-','color',cmap(i,:))


    h = headGains(ii,end);
    pp = -PIscores(ii,end);
    [pfit, S] = polyfit(h, pp, 1);
    R_squared = 1 - (S.normr/norm(pp - mean(pp)))^2;

    subplot(2,3,4);hold on;
    plot(R_squared,i,'o','color',cmap(i,:),'markersize',12);
    ylim([0,6])
    xlim([0,.25])
    set(gca,'fontsize',16)
    xlabel('R^2')
    ylabel('goal strength')

    if i==5
        subplot(2,3,3);hold on;
        
        set(gca,'fontsize',16)
        plot(h,pp,'o','markeredgecolor','none','markerfacecolor','k')
        
        xx = linspace(0,.6,100);
        plot(xx,pfit(1).*xx+pfit(2),'--k');
        xlabel('final map strength')
        ylabel('final preference')

        subplot(2,3,5);hold on;
        plot(mean(-PIscores(ii,:)),'-k');
        set(gca,'fontsize',16)
        xlabel('time (a.u.)')
        ylabel('strength of behavioral preference')

    end
end

if plotMaps
    nx = p.N;
    thetaF  = linspace(0,2*pi,nx+1);
    theta  = thetaF(1:end-1);
    dtheta = thetaF(2)-thetaF(1);

    [~,~,~,~,~,~,~,~,~,res] = runLearningSims(p,.5,17,params);
    nx = p.N;
    headingLabels = 1:nx;
    arenaLabels   = fliplr(headingLabels);
    iter = [1,50,min(300,size(res.weights.heading,3))];

    figure;hold on;set(gcf,'Position',[200 200 1000 600],'color','w')
    for i=1:numel(iter)
        hMat = squeeze(res.weights.heading(:,:,iter(i)));
        subplot(2,numel(iter),i);imagesc(theta+dtheta/2,theta+dtheta/2,hMat);colormap(gray);
        caxis([0,1]);hold on;
        plot([pi,pi],[0,2*pi],'--k')
        set(gca,'fontsize',16,'ydir','normal')
        xlabel('compass heading')
        ylabel('arena heading ')

        pJ = [];
        for k = 1:nx
            h = arenaLabels(k);
            candArenaHeadings = [circIndex(h-1:h+1,nx),circIndex(h+nx/2-1:h+nx/2+1,nx)];
            activities = [.5,1,.5,.5,1,.5];
            dp1 = activities*hMat(candArenaHeadings,k)./numel(activities);
            dp2 = activities*hMat(candArenaHeadings,circIndex(k-nx/2,nx))./numel(activities);
            xx = [1./dp1,1./dp2];
            xx = xx./sum(xx);
            pJ(k) = xx(1); %this is pStay
        end

        rH    = sum(pJ.*exp(1i.*theta));
        angH(i)  = angle(rH);
        goalH(i) = mod(angH(i),2*pi);
        strengthH(i) = abs(rH)./12;
        
        
        subplot(2,numel(iter),i+numel(iter));hold on;
        plot(1-pJ);
        plot([1,32],[.5,.5],'--k');
        plot([17,17],[0,1],'--b');
        xlim([1,32]);
        ylim([.2,.8]);
        set(gca,'fontsize',16)
        xlabel('compass heading (wedges)')
        ylabel('p(jump)')
        
    end
end
