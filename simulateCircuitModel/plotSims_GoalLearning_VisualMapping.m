function plotSims_GoalLearning_VisualMapping(plotNum,headLocs0,headGains0,goalLocs0,goalGains0,distDanger0,PIscores0,angdiff0,gain0)
% PLOTSIMS_GOALLEARNING_VISUALMAPPING plots the output of simulations
% designed to probe the evolution of the visual mapping and the goal
% heading
%
% INPUTS:
%   plotting key for plotNum
%       1: vector maps
%       2: individual trajectories
%       3: misalignment across individual flies
%       4: variance in goal updates
%       5: fraction of flies with different initial goal headings
%       6: average learning trajectories, split based on initial strength of goal heading
%       7: CDF of initial conditions
%   headLocs0:      location of most stable compass heading over time
%   headGains0:     strength of most stable compass heading over time
%   goalLocs0:      location of goal heading over time
%   goalGains0:     strength of goal heading over time
%   distDanger0:    angular distance between goal heading and center of safe 
%                       zone over time 
%   PIscores0:      PI score over time
%   angDiff0:       angular distance between most stable compass heading and
%                       goal heading
%   gain0:          array containing initial strengths of goal heading
%

dI = 200;
I0 = 0;
gain0 = reshape(gain0,[dI*dI,5]);
gain0 = gain0(1:dI,:);
gain0 = gain0(:);
distDanger = distDanger0(I0+1:end,:);
goalGains  = goalGains0(I0+1:end,:);
goalLocs   = goalLocs0(I0+1:end,:);
headGains  = headGains0(I0+1:end,:);
headLocs   = headLocs0(I0+1:end,:);
PIscores   = PIscores0(I0+1:end,:);
angdiff    = angdiff0(I0+1:end,:);
gain       = gain0(I0+1:end,:);

g = unique(gain0);
[niter,nt] = size(goalGains);

%compute variance in goal updates
g0 = goalLocs(:,1:end-1);
gf = repmat(goalLocs(:,end),[1,nt-1]);
g1 = gf+2*pi;
g2 = gf;
g3 = gf-2*pi;
dG  = cat(3,g1,g2,g3)-cat(3,g0,g0,g0);
[~,inds] = min(abs(dG),[],3);
for i=1:niter
    for j=1:nt-1
        min_dG(i,j) = dG(i,j,inds(i,j));
    end
end

h0 = min_dG(:,1:end-1);
h1 = min_dG(:,2:end)-2*pi;
h2 = min_dG(:,2:end);
h3 = min_dG(:,2:end)+2*pi;
vG  = cat(3,h1,h2,h3)-cat(3,h0,h0,h0);
[~,inds] = min(abs(vG),[],3);
for i=1:niter
    for j=1:nt-2
        min_vG(i,j) = vG(i,j,inds(i,j));
    end
end

G = nan(niter,800);
A = nan(niter,800);
P = nan(niter,800);

dG = nan(niter,800);
vG = nan(niter,800);


cmapA = viridis(100);
cmapC = flipud(inferno(105));
dc = (cmapC(1,:)-cmapA(end,:))./5;
cmapB = cmapA(end,:)+dc;
for i=2:4
    cmapB = [cmapB;cmapA(end,:)+i*dc];
end
redblueI = [cmapA;cmapB;cmapC];
for i=1:niter
    [~,ii] = min((goalGains(i,:)));
    G(i,400-ii+1:400-ii+nt) = goalGains(i,:);
    A(i,400-ii+1:400-ii+nt) = angdiff(i,:);
    P(i,400-ii+1:400-ii+nt) = PIscores(i,:);%-PIscores(i,1);
    
    dG(i,400-ii+1:400-ii+nt-1) = min_dG(i,:);
    vG(i,400-ii+1:400-ii+nt-2) = min_vG(i,:);
end

%sort by depth of gain drop over short window :
[~,jj] = sort(nanmean(G(:,401:410),2));
imin = 1:dI:5*dI;


%------------------------ generate vector maps ---------------------------%
if ismember(1,plotNum)

    cmapA = viridis(100);
    cmapC = flipud(inferno(105));
    dc = (cmapC(1,:)-cmapA(end,:))./5;
    cmapB = cmapA(end,:)+dc;
    for i=2:4
        cmapB = [cmapB;cmapA(end,:)+i*dc];
    end
    cmap = [cmapA;cmapB;cmapC];
    
    cmap = flipud(inferno(100));
    nc = size(cmap,1);
    
    
    headLocs = headLocs+pi;
    headLocs(headLocs>2*pi) = headLocs(headLocs>2*pi)-2*pi;
    goalLocs = goalLocs+pi;
    goalLocs(goalLocs>2*pi) = goalLocs(goalLocs>2*pi)-2*pi;
    
    hL1 = headLocs(201:end,1:end-1);
    hL2 = headLocs(201:end,2:end);
    gL1 = goalLocs(201:end,1:end-1);
    gL2 = goalLocs(201:end,2:end);
    
    h1 = headGains(201:end,1:end-1);
    h2 = headGains(201:end,2:end);
    g1 = goalGains(201:end,1:end-1);
    g2 = goalGains(201:end,2:end);
    d1 = distDanger(201:end,1:end-1);
    d2 = distDanger(201:end,2:end);
    p1 = PIscores(201:end,1:end-1);
    
    gL1 = gL1(:);
    gL2 = gL2(:);
    hL1 = hL1(:);
    hL2 = hL2(:);
    
    d1 = d1(:);
    d2 = d2(:);
    g1 = g1(:);
    g2 = g2(:);
    h1 = h1(:);
    h2 = h2(:);
    p1 = p1(:);
    
    nbins = 15;
    dbins = linspace(0,pi/2,nbins);
    gbins = linspace(.1,.8,nbins);
    hbins = linspace(.1,.6,nbins);
    
    nmin = 50;
    
    nall = [];
    figure;set(gcf,'Position',[200 200 1400 400],'color','w');hold on;
    subplot(1,3,3);hold on;
    
    for i=1:nbins-1
        for j=1:nbins-1
            ii = find(g1>gbins(j) & g1<gbins(j+1) & d1>dbins(i) & d1<dbins(i+1));
            nsamp(i,j) = numel(ii);
            pitmp(i,j) = mean(p1(ii));
        end
    end
    nall = cat(3,nall,nsamp);
    pitmp(nsamp<nmin) = nan;
    cvec  = linspace(nanmin(pitmp(:)),nanmax(pitmp(:)),nc);
    for i=1:nbins-1
        for j=1:nbins-1
            ii = find(g1>gbins(j) & g1<gbins(j+1) & d1>dbins(i) & d1<dbins(i+1));
            [~,jj] = min(abs(mean(p1(ii))-cvec));
            if nsamp(i,j)>nmin
                dist = sqrt(mean(g2(ii)-g1(ii)).^2 + mean(d2(ii)-d1(ii)).^2);
                dist = 10*dist.^.75;
                quiver(mean(g1(ii)),mean(d1(ii)),mean(g2(ii)-g1(ii))./dist,mean(d2(ii)-d1(ii))./dist,2,'color',cmap(jj,:),'linewidth',2,'maxheadsize',10);
            end
        end
    end
    set(gca,'fontsize',16)
    xlim([.1,.8]);
    ylim([0,pi/2])
    yticks([0,pi/4,pi/2])
    title(['colormap (PI score): ',num2str(min(pitmp(:))),', ',num2str(max(pitmp(:)))])
    yticklabels({'0','\pi/4','\pi/2'})
    ylabel('distance to safety')
    xlabel('goal gain')
    
    
    subplot(1,3,2);hold on;
    for i=1:nbins-1
        for j=1:nbins-1
            ii = find(g1>gbins(j) & g1<gbins(j+1) & h1>hbins(i) & h1<hbins(i+1));
            nsamp(i,j) = numel(ii);
            pitmp(i,j) = mean(p1(ii));
        end
    end
    nall = cat(3,nall,nsamp);
    pitmp(nsamp<nmin) = nan;
    cvec  = linspace(nanmin(pitmp(:)),nanmax(pitmp(:)),nc);
    for i=1:nbins-1
        for j=1:nbins-1
            ii = find(g1>gbins(j) & g1<gbins(j+1) & h1>hbins(i) & h1<hbins(i+1));
            [~,jj] = min(abs(mean(p1(ii))-cvec));        
            if nsamp(i,j)>nmin
                dist = sqrt(mean(g2(ii)-g1(ii)).^2 + mean(h2(ii)-h1(ii)).^2);
                dist =  15.*dist.^.75;%15.*dist.^.6;
                quiver(mean(g1(ii)),mean(h1(ii)),mean(g2(ii)-g1(ii))./dist,mean(h2(ii)-h1(ii))./dist,2,'color',cmap(jj,:),'linewidth',2,'maxheadsize',10);
            end
        end
    end
    plot([0.1,0.8],[0.1,0.8],'--k')
    xlim([.1,.8])
    ylim([.1,.6]);
    xticks(.1:.1:.8)
    yticks(.1:.1:.6)
    title(['colormap (PI score): ',num2str(min(pitmp(:))),', ',num2str(max(pitmp(:)))])
    set(gca,'fontsize',16)
    ylabel('strength of most stable compass heading')
    xlabel('strength of goal heading')
    
    
    hLbins = linspace(0,2*pi,nbins);
    gLbins = linspace(0,2*pi,nbins);
    subplot(1,3,1);hold on;
    for i=1:nbins-1
        for j=1:nbins-1
            ii = find(gL1>gLbins(j) & gL1<gLbins(j+1) & hL1>hLbins(i) & hL1<hLbins(i+1));
            nsamp(i,j) = numel(ii);
            pitmp(i,j) = mean(p1(ii));
        end
    end
    nall = cat(3,nall,nsamp);
    pitmp(nsamp<nmin) = nan;
    cvec  = linspace(nanmin(pitmp(:)),nanmax(pitmp(:)),nc);
    for i=1:nbins-1
        for j=1:nbins-1
            ii = find(gL1>gLbins(j) & gL1<gLbins(j+1) & hL1>hLbins(i) & hL1<hLbins(i+1));
            [~,jj] = min(abs(mean(p1(ii))-cvec));
            if nsamp(i,j)>nmin
                dist = sqrt(mean(gL2(ii)-gL1(ii)).^2 + mean(hL2(ii)-hL1(ii)).^2);
                dist = 5.*dist.^.75;
                quiver(mean(gL1(ii)),mean(hL1(ii)),mean(gL2(ii)-gL1(ii))./dist,mean(hL2(ii)-hL1(ii))./dist,3,'color',cmap(jj,:),'linewidth',2,'maxheadsize',10);
                caxis([-11,-2]);
            end
        end
    end
    plot([0,2*pi],[0,2*pi],'--k')
    xlim([0,2*pi])
    ylim([0,2*pi]);
    xticks([0,pi,2*pi])
    yticks([0,pi,2*pi])
    xticklabels({'0','\pi','2\pi'})
    yticklabels({'0','\pi','2\pi'})
    title(['colormap (PI score): ',num2str(min(pitmp(:))),', ',num2str(max(pitmp(:)))])
    set(gca,'fontsize',16)
    ylabel('most stable compass heading')
    xlabel('goal heading')
    
    
    figure;hold on;set(gcf,'Position',[200 200 1400 400],'color','w')
    xlabels = {'goal strength', 'goal orientation', 'goal strength'};
    ylabels = {'map strength', 'map orientation', 'distance to safety'};
    for i=1:3
        Poccur = log(nall(:,:,i)./(size(headGains,1).*size(headGains,2)));
        subplot(1,3,i);imagesc(Poccur);set(gca,'ydir','normal');colormap(flipud(gray));
        caxis([-11,-2]);
        set(gca,'fontsize',16)
        xlabel(xlabels{i})
        ylabel(ylabels{i})
        title(['colormap (PI score): ',num2str(min(Poccur(:))),', ',num2str(max(Poccur(:)))])
    end


end


%---------------- plot individual trajectories ---------------------------%

if ismember(2,plotNum) 
    figure;set(gcf,'Position',[200 200 1400 400],'color','w')
    xtx = 0:.2:.8;
    ytx = [0,pi/4,pi/2];
    for i=1:size(PIscores,1)
        tf1   = find(goalGains(i,:)<0.6,1,'last');
        tf2   = find(distDanger(i,:)>pi/8,1,'last');
        if numel(tf1)>0 && numel(tf2)>0
            tf(i) = max(tf1,tf2)./size(PIscores,2);
        else
            tf(i) = nan;
        end
    end
    
    for i=1:5
        inds = 1:2:size(PIscores,2);
        ii = find(gain==g(i));
    
        [~,jj] = sort(mean(PIscores(ii,:),2));
    
        subplot(2,5,i);hold on;
        k = randperm(10,1);
        plot(goalGains(ii(jj(k)),inds),distDanger(ii(jj(k)),inds),'-k')
        scatter(goalGains(ii(jj(k)),inds),distDanger(ii(jj(k)),inds),100,pi-angdiff(ii(jj(k)),inds),'filled')
        caxis([0,pi]);xlim([0,.8]);
        title(['niter = ',num2str(round(tf(ii(jj(k))),2))])
        colormap(flipud(inferno))
        yticks(ytx)
        xticks(xtx)
        xticklabels({})
        yticklabels({})
        set(gca,'fontsize',16)
    
        [~,jj] = sort(mean(PIscores(ii,:),2),'descend');
        subplot(2,5,5+i);hold on;
        k = randperm(10,1);
        plot(goalGains(ii(jj(k)),inds),distDanger(ii(jj(k)),inds),'-k')
        scatter(goalGains(ii(jj(k)),inds),distDanger(ii(jj(k)),inds),100,pi-angdiff(ii(jj(k)),inds),'filled')
        caxis([0,pi]);xlim([0,.8])
        title(['niter = ',num2str(round(tf(ii(jj(k))),2))])
        colormap(flipud(inferno))
        yticks(ytx)
        xticks(xtx)
        xticklabels({})
        yticklabels({})
        set(gca,'fontsize',16)
    end
        
    for i=1:1000
        for j=1:322
            distDangerHead0(i,j) = min(abs(headLocs0(i,j) - [pi/2,3*pi/2]));
        end
    end
end


%----------------plot misalignment across individual flies----------------%
if ismember(3,plotNum) 
    figure;hold on;set(gcf,'Position',[200 200 1000 800],'color','w')
    set(gcf,'color','w')
    
    imin = 1:dI:5*dI;
    cmap = inferno(6);
    cmap = flipud(cmap(1:end-1,:));
    for i=1:numel(imin)
        ii = imin(i):imin(i)+dI-1;
    
        subplot(2,3,4);hold on;
        plot(nanmean(G(jj(ii),351:end)),'linewidth',2,'color',cmap(i,:));
        if i==numel(imin)
            plot([50,50],[0,.8],'--k','linewidth',1);
        end
        xlim([0,300]);
        ylim([0,.8])
        set(gca,'fontsize',16);
        ylabel('goal gain');
        xlabel('time')

        subplot(2,3,5);hold on;
        plot(nanmean(A(jj(ii),351:end)),'linewidth',2,'color',cmap(i,:));
        if i==numel(imin)
            plot([50,50],[0,pi],'--k','linewidth',1);
        end
    
        xlim([0,300]);
        set(gca,'fontsize',16);
        ylabel('misalignment');
        xlabel('time')

        subplot(2,3,6);hold on;
        plot(nanmean(P(jj(ii),351:end)),'linewidth',2,'color',cmap(i,:));
        if i==numel(imin)
            plot([50,50],[-.3,.6],'--k','linewidth',1);
        end
        xlim([0,300]);
        ylim([-.3,.6])
        set(gca,'fontsize',16);
        ylabel('PI scores');
        xlabel('time')
        
        c = histcounts(gain(jj(ii)),.05:.1:.55);
        c(c<1) = .01;
        c = c./numel(ii);
        labels = {'','','','',''};
        
    end
    
    ax2 = subplot(2,3,1);hold on;
    imagesc(G(jj,351:end),'AlphaData',~isnan(G(jj,351:end)));
    caxis([0,.8])
    colormap(ax2,flipud(inferno))
    plot([50,50],[0,500],'--k','linewidth',1);
    xlim([0,300]);
    ylim([0,500]);
    set(gca,'fontsize',16);
    yticklabels({})
    title('goal gain');
    xlabel('time')
    
    ax3 = subplot(2,3,2);hold on;
    imagesc(pi-A(jj,351:end),'AlphaData',~isnan(A(jj,351:end)));
    colormap(ax3,flipud(inferno))
    plot([50,50],[0,500],'--k','linewidth',1);
    xlim([0,300]);
    ylim([0,500]);
    set(gca,'fontsize',16);
    yticklabels({})
    title('misaslignment');
    xlabel('time')
    
    ax4 = subplot(2,3,3);hold on;
    imagesc(P(jj,351:end),'AlphaData',~isnan(P(jj,351:end)));
    caxis([-.7,.7])
    colormap(ax4,redblueI)
    plot([50,50],[0,500],'--k','linewidth',1);
    xlim([0,300]);
    ylim([0,500]);
    set(gca,'fontsize',16);
    yticklabels({})
    title('PI scores');
    xlabel('time')
end



%--------------------- plot variance in goal updates ---------------------%
if ismember(4,plotNum) 
    figure;hold on;set(gcf,'Position',[200 200 1000 800],'color','w')
    set(gcf,'color','w')
    
    vvG = (movvar(vG,10,0,2));
    
    cmap = inferno(6);
    cmap = flipud(cmap(1:end-1,:));
    for i=1:numel(imin)
        ii = imin(i):imin(i)+dI-1;
    
        subplot(2,3,4);hold on;
        plot(nanmean(abs(dG(jj(ii),351:end))),'linewidth',2,'color',cmap(i,:));
        if i==numel(imin)
            plot([50,50],[0,2],'--k','linewidth',1);
        end
        xlim([0,300]);
        %ylim([0,.8])
        set(gca,'fontsize',16);
        ylabel('distance from final goal');
    
        subplot(2,3,5);hold on;
        plot(nanmean(log10(vvG(jj(ii),351:end))),'linewidth',2,'color',cmap(i,:));
        if i==numel(imin)
            plot([50,50],[-8,-1],'--k','linewidth',1);
        end
    
        xlim([0,300]);
        set(gca,'fontsize',16);
        ylabel('change in distance');
        
        subplot(2,3,6);hold on;
        errorbar(i,nanmean(log10(nanvar(vG(jj(ii),351:end),[],2))),nanstd(log(nanvar(vG(jj(ii),351:end),[],2))),...
            'o','linewidth',2,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'markeredgecolor','none');
        xlim([0,6])
        set(gca,'fontsize',16);
        ylabel('variance in goal updates');
        
    end
    
    
    ax2 = subplot(2,3,1);hold on;
    imagesc(abs(dG(jj,351:end)),'AlphaData',~isnan(dG(jj,351:end)));
    colormap(ax2,inferno)
    plot([50,50],[0,500],'--k','linewidth',1);
    xlim([0,300]);
    ylim([0,500]);
    set(gca,'fontsize',16);
    yticklabels({})
    title('distance from final goal');
    xlabel('time')
    
    ax3 = subplot(2,3,2);hold on;
    imagesc(log10(vvG(jj,351:end)),'AlphaData',~isnan(vvG(jj,351:end)));
    caxis([-3,max(log10(vvG(:)))]);
    colormap(ax3,inferno)
    plot([50,50],[0,500],'--k','linewidth',1);
    xlim([0,300]);
    ylim([0,500]);
    set(gca,'fontsize',16);
    yticklabels({})
    title('change in distance');
    xlabel('time')
    
    ax4 = subplot(2,3,3);hold on;
    plot(log10(nanvar(vG(jj,:),[],2)),1:niter,'-k')
    set(gca,'fontsize',16);
    yticklabels({})
    title('variance in goal updates');
    xlabel('time')
end


%--------------plot fraction of model flies w/ different gains------------%
if ismember(5,plotNum) 
    figure;hold on;set(gcf,'Position',[200 200 1400 300],'color','w')
    cmap = gray(6);
    cmap = flipud(cmap(1:end-1,:));
    for i=1:numel(imin)
        ii = imin(i):imin(i)+dI-1;
    
        c = histcounts(gain(jj(ii)),.05:.1:.55);
        c(c<1) = .01;
        c = c./numel(ii);
        labels = {'','','','',''};
        
        ax1 = subplot(1,5,i);hold on;
        pie(c,labels);
        fill([-1,1,1,-1,-1],[1.2,1.2,1.5,1.5,1.3],cmap(i,:))
        axis off
        colormap(ax1,cmap)
        set(gca,'fontsize',16)
        if i==3
            title('--- minimum gain -->')
        end
        if i==5
            legend({'initial gain = 0.1','initial gain = 0.2',...
                'initial gain = 0.3','initial gain = 0.4','initial gain = 0.5'})
        end
    end
end


if ismember(6,plotNum)

    distDanger = distDanger*180./pi;
    nt = size(PIscores,2);
    nperms = 20;
    
    figure; hold on;set(gcf,'Position',[200 200 400 800],'color','w')
    % ----------------------------------------------------------------------- %
    [~,jj] = sort(goalGains(:,1));
    
    PIs1 = [];PIs2 = [];PIm1 = [];PIm2 = [];
    DDs1 = [];DDs2 = [];DDm1 = [];DDm2 = [];
    Gs1 = []; Gs2 = []; Gm1 = []; Gm2 = [];
    for i=1:100
        ii = randperm(1000,2*nperms);
        i1 = find(goalGains(jj(ii),1)< nanmedian(goalGains(jj(ii),1)));
        i2 = find(goalGains(jj(ii),1)>=nanmedian(goalGains(jj(ii),1)));
    
        PIs1 = [PIs1;std(PIscores(jj(ii(i1)),:))./sqrt(nperms)];
        PIm1 = [PIm1;mean(PIscores(jj(ii(i1)),:))];
        PIs2 = [PIs2;std(PIscores(jj(ii(i2)),:))./sqrt(nperms)];
        PIm2 = [PIm2;mean(PIscores(jj(ii(i2)),:))];
        
        Gs1 = [Gs1;std(goalGains(jj(ii(i1)),:))./sqrt(nperms)];
        Gm1 = [Gm1;mean(goalGains(jj(ii(i1)),:))];
        Gs2 = [Gs2;std(goalGains(jj(ii(i2)),:))./sqrt(nperms)];
        Gm2 = [Gm2;mean(goalGains(jj(ii(i2)),:))];
        
        DDs1 = [DDs1;std(distDanger(jj(ii(i1)),:))./sqrt(nperms)];
        DDm1 = [DDm1;mean(distDanger(jj(ii(i1)),:))];
        DDs2 = [DDs2;std(distDanger(jj(ii(i2)),:))./sqrt(nperms)];
        DDm2 = [DDm2;mean(distDanger(jj(ii(i2)),:))];
    end
    
    
    m1 = mean(PIm1);m2 = mean(PIm2);
    s1 = std(PIm1);s2 = std(PIm2);
    subplot(3,1,3);hold on;
    fill([1:nt,nt:-1:1,1],[m1+s1,fliplr(m1-s1),m1(1)+s1(1)],[.8,.8,.8],'edgecolor','none')
    plot(m1);
    fill([1:nt,nt:-1:1,1],[m2+s2,fliplr(m2-s2),m2(1)+s2(1)],[.6,.6,.6],'edgecolor','none')
    plot(m2);
    xlim([1,300])
    ylim([-.1,.5])
    set(gca,'fontsize',16)
    ylabel('time (a.u.)')
    xlabel('PI score')
    
    m1 = mean(DDm1);m2 = mean(DDm2);
    s1 = std(DDm1);s2 = std(DDm2);
    subplot(3,1,2);hold on;
    fill([1:nt,nt:-1:1,1],[m1+s1,fliplr(m1-s1),m1(1)+s1(1)],[.8,.8,.8],'edgecolor','none')
    plot(m1);
    fill([1:nt,nt:-1:1,1],[m2+s2,fliplr(m2-s2),m2(1)+s2(1)],[.6,.6,.6],'edgecolor','none')
    plot(m2);
    xlim([1,300])
    ylim([0,52.5])
    yticks([0,15,30,45])
    set(gca,'fontsize',16)
    ylabel('time (a.u.)')
    xlabel('distance to safety')

    % ----------------------------------------------------------------------- %
    
    m1 = mean(Gm1);m2 = mean(Gm2);
    s1 = std(Gm1);s2 = std(Gm2);
    
    subplot(3,1,1);hold on;
    fill([1:nt,nt:-1:1,1],[m1+s1,fliplr(m1-s1),m1(1)+s1(1)],[.8,.8,.8],'edgecolor','none')
    plot(m1);
    fill([1:nt,nt:-1:1,1],[m2+s2,fliplr(m2-s2),m2(1)+s2(1)],[.6,.6,.6],'edgecolor','none')
    plot(m2);
    xlim([1,300])
    ylim([.1,.8])
    set(gca,'fontsize',16)
    ylabel('time (a.u.)')
    xlabel('goal strength')
    legend({'','weak initial goal','','strong initial goal'})

end

if ismember(7,plotNum)
    figure; hold on;set(gcf,'Position',[200 200 400 800],'color','w')
    subplot(2,1,1);plot(sort(goalGains(:,1)),linspace(0,1,size(goalGains,1)),'-k','linewidth',2)
    xlim([-.2,1])
    ylim([0,1])
    xlabel('goal strength')
    ylabel('CDF')
    set(gca,'fontsize',16)

    subplot(2,1,2);plot(sort(goalLocs(:,1)),linspace(0,1,size(goalLocs,1)),'-k','linewidth',2)
    xlim([0,2*pi])
    ylim([0,1])
    xlabel('goal location')
    ylabel('CDF')
    set(gca,'fontsize',16)

end

end





