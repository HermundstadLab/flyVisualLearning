function pvals = plotLearningEvolution(plotNum,sL,rsL,sNL,rsNL,sKC,rsKC,p,combineData)
% PLOTLEARNINGEVOLUTION plots different features of learning data
%
% INPUTS:
%   plotNum: integer value that specifies which type of plot to generate:
%       1: plot distance to safety for laser-trained and no-laser control flies
%       2: plot learning evolution, split into 2 groups
%       3: plot learning evolution, split into 3 groups
%   rsL, rsNL, rsKC: data structures containing the raw heading trajectories
%   sL, sNL, sKC:  data structures containining segmented data
%   p: structure containing parameter values
%   combineData: logical that specifies whether or not to combine datasets.
%       Default is false.
%


if nargin<9
    combineData = false;
end


nx = p.nx/2;

%--------------------- plot distance to safety ---------------------------%
if ismember(1,plotNum)
    ntrials = numel(sL.fly(1).trial);
    for i=1:numel(sL.fly)
        for j=1:ntrials
            ddanger(i,j)   = min(abs(sL.fly(i).trial(j).shiftedFull.prefLoc-[nx/2-1,nx/2,3*nx/2-1,3*nx/2]));
        end
    end
    
    for i=1:numel(sNL.fly)
        for j=1:ntrials
            ddangerNL(i,j) = min(abs(sNL.fly(i).trial(j).shiftedFull.prefLoc-[nx/2-1,nx/2,3*nx/2-1,3*nx/2]));
        end
    end
    
    ddanger   = ddanger.*360./p.nx; 
    ddangerNL = ddangerNL.*360./p.nx; 
    
    for i=1:ntrials
        pvals(i) = ranksum(ddanger(:,i),ddangerNL(:,i));
    end
    
    figure;set(gcf,'Position',[200 200 600 400],'color','w');hold on;
    indsL  = (1:ntrials)-.1;
    indsNL = (1:ntrials)+.1;
    errorbar(indsL, nanmean(ddanger,  1),nanstd(ddanger,  [],1)./sqrt(sum(1-isnan(ddanger)  )),'ko-','LineWidth',2,'markersize',12)
    errorbar(indsNL,nanmean(ddangerNL,1),nanstd(ddangerNL,[],1)./sqrt(sum(1-isnan(ddangerNL))),'o-','LineWidth',2,'markersize',12,'color',[.6,.6,.6])
    set(gca,'fontsize',16,'tickdir','out')
    ylabel('distance to safety (deg)')
    xlabel('trial')
    legend({'laser-trained','no-laser controls'})
    xticks(1:9)
    xlim([0,10])
end

%--------------- plot learning evolution, split into 2 groups ------------%
if ismember(2,plotNum)
    
    if combineData
        s = {sL,sKC};
        rs = {rsL,rsKC};
    else
        s = {sL};
        rs = {rsL};
    end
    

    fly = 1;
    for k=1:numel(s)
        stmp = s{k};
        rstmp = rs{k};
        for i=1:numel(stmp.fly)
            for j=1:numel(stmp.fly(1).trial)
                x  = rstmp.fly(i).trial(j).x;
                xw = wrap(x,rstmp.fly(i).punishedBar);
                
                inds = stmp.fly(i).trial(j).PIinds; 
                xw = xw(inds);
                
                ddanger(fly,j)   = min(abs(stmp.fly(i).trial(j).shiftedFull.prefLoc-[23,24,71,72]));      
                goalLoc(fly,j) = stmp.fly(i).trial(j).shiftedFull.prefLoc;       
                x  = quickWrap(xw + stmp.fly(i).trial(j).shifted.shift,nx);
                r(fly,j) = (numel(find(x>=12 & x<36))-numel([find(x<12);find(x>=36)]))./numel(x);
                PI(fly,j)    = stmp.fly(i).trial(j).PI;
            end
            fly = fly+1;
        end
    end
    
    ii = find(nanmean(r(:,2),2)<nanmedian(r(:,2)));
    jj = find(nanmean(r(:,2),2)>=nanmedian(r(:,2)));
    
    pvals = nan(1,p.nTrials);
    for i=2:p.nTrials
        pvals(1,i) = ranksum(r(jj,i),r(ii,i));
        pvals(2,i) = ranksum(ddanger(jj,i),ddanger(ii,i));
        pvals(3,i) = ranksum(PI(jj,i),PI(ii,i));
    end
    
    figure;set(gcf,'Position',[200 200 1200 400],'color','w');hold on;hold on
    i0 = 2;
    inds1 = 2:p.nTrials;
    inds2 = inds1+.2;
    inds3 = inds2+.2;
    
    cS = [0,0,0];
    cW = [252,165,10]./255;
    
    subplot(1,3,1);hold on;
    errorbar(inds1,nanmean(r(ii,i0:end),1),nanstd(r(ii,i0:end),[],1)./sqrt(sum(1-isnan(r(ii,i0:end)))),'color',cW,'marker','none','linestyle','none','LineWidth',0.75)
    errorbar(inds2,nanmean(r(jj,i0:end),1),nanstd(r(jj,i0:end),[],1)./sqrt(sum(1-isnan(r(jj,i0:end)))),'color',cS,'marker','none','linestyle','none','linewidth',0.75)
    plot(inds1,nanmean(r(ii,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cW,'color',cW,'LineWidth',0.75)
    plot(inds2,nanmean(r(jj,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cS,'color',cS,'LineWidth',0.75)
    ylabel('gain of behavior')
    xlabel('trial')
    set(gca,'fontsize',16)
    xticks(inds1)
    xlim([inds1(1)-1,inds1(end)+1])
    ylim([.1,.9]);
    yticks(.1:.1:.9);
    
    subplot(1,3,2);hold on;
    errorbar(inds1,nanmean(ddanger(ii,i0:end),1),nanstd(ddanger(ii,i0:end),[],1)./sqrt(sum(1-isnan(ddanger(ii,i0:end)))),'color',cW,'marker','none','linestyle','none','LineWidth',0.75)
    errorbar(inds2,nanmean(ddanger(jj,i0:end),1),nanstd(ddanger(jj,i0:end),[],1)./sqrt(sum(1-isnan(ddanger(jj,i0:end)))),'color',cS,'marker','none','linestyle','none','linewidth',0.75)
    plot(inds1,nanmean(ddanger(ii,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cW,'color',cW,'LineWidth',0.75)
    plot(inds2,nanmean(ddanger(jj,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cS,'color',cS,'LineWidth',0.75)
    ylabel('distance to safety')
    xlabel('trial')
    set(gca,'fontsize',16)
    xticks(inds1)
    xlim([inds1(1)-1,inds1(end)+1])
    ylim([0,16]);
    yticks([0,4,8,12,16]);
    yticklabels({num2str(0),num2str(4*90./24),num2str(8*90./24),num2str(12*90./24),num2str(16*90./24)});
    
    subplot(1,3,3);hold on;
    errorbar(inds1,nanmean(PI(ii,i0:end),1),nanstd(PI(ii,i0:end),[],1)./sqrt(sum(1-isnan(PI(ii,i0:end)))),'color',cW,'marker','none','linestyle','none','LineWidth',0.75)
    errorbar(inds2,nanmean(PI(jj,i0:end),1),nanstd(PI(jj,i0:end),[],1)./sqrt(sum(1-isnan(PI(jj,i0:end)))),'color',cS,'marker','none','linestyle','none','linewidth',0.75)
    plot(inds1,nanmean(PI(ii,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cW,'color',cW,'LineWidth',0.75)
    plot(inds2,nanmean(PI(jj,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cS,'color',cS,'LineWidth',0.75)
    ylabel('PI score')
    xlabel('trial')
    set(gca,'fontsize',16)
    xticks(inds1)
    xlim([inds1(1)-1,inds1(end)+1])
    legend({'','','weak initial goal','strong initial goal'});

end


%--------------- plot learning evolution, split into 3 groups ------------%
if ismember(3,plotNum)
    
    if combineData
        s = {sL,sKC};
        rs = {rsL,rsKC};
    else
        s = {sL};
        rs = {rsL};
    end
    
    fly = 1;
    for k=1:numel(s)
        stmp = s{k};
        rstmp = rs{k};
        for i=1:numel(stmp.fly)
            for j=1:numel(stmp.fly(1).trial)
                x  = rstmp.fly(i).trial(j).x;
                xw = wrap(x,rstmp.fly(i).punishedBar);
                
                inds = stmp.fly(i).trial(j).PIinds; 
                xw = xw(inds);
                
                ddanger(fly,j)   = min(abs(stmp.fly(i).trial(j).shiftedFull.prefLoc-[23,24,71,72]));
                
                goalLoc(fly,j) = stmp.fly(i).trial(j).shiftedFull.prefLoc;
                
                x  = quickWrap(xw + stmp.fly(i).trial(j).shifted.shift,nx);
                r(fly,j) = (numel(find(x>=12 & x<36))-numel([find(x<12);find(x>=36)]))./numel(x);
                PI(fly,j)    = stmp.fly(i).trial(j).PI;
            end
            fly = fly+1;
        end
    end
    
    %to consider three groups
    q = quantile(r(:,2),[0,1/3,2/3,1]);
    ii = find(r(:,2)<q(2));
    jj = find(r(:,2)>=q(2) & r(:,2)<q(3));
    kk = find(r(:,2)>=q(3));
    
    
    pvals    = nan(1,p.nTrials);
    for i=2:p.nTrials
        pvals(1,i) = ranksum(r(kk,i),r(ii,i));
        pvals(2,i) = ranksum(ddanger(kk,i),ddanger(ii,i));
        pvals(3,i) = ranksum(PI(kk,i),PI(ii,i));
    end
    
    figure;set(gcf,'Position',[200 200 1200 400],'color','w');hold on;hold on
    i0 = 2;
    inds1 = 2:9;
    inds2 = inds1+.2;
    inds3 = inds2+.2;
    
    cS = [0,0,0];
    cM = [147,38,103]./255;
    cW = [252,165,10]./255;

    subplot(1,3,1);hold on;
    errorbar(inds1,nanmean(r(ii,i0:end),1),nanstd(r(ii,i0:end),[],1)./sqrt(sum(1-isnan(r(ii,i0:end)))),'color',cW,'marker','none','linestyle','none','LineWidth',0.75)
    errorbar(inds2,nanmean(r(jj,i0:end),1),nanstd(r(jj,i0:end),[],1)./sqrt(sum(1-isnan(r(jj,i0:end)))),'color',cM,'marker','none','linestyle','none','linewidth',0.75)
    errorbar(inds3,nanmean(r(kk,i0:end),1),nanstd(r(kk,i0:end),[],1)./sqrt(sum(1-isnan(r(kk,i0:end)))),'color',cS,'marker','none','linestyle','none','linewidth',0.75)
    plot(inds1,nanmean(r(ii,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cW,'color',cW,'LineWidth',0.75)
    plot(inds2,nanmean(r(jj,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cM,'color',cM,'LineWidth',0.75)
    plot(inds3,nanmean(r(kk,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cS,'color',cS,'LineWidth',0.75)
    ylabel('gain of behavior')
    xlabel('trial')
    set(gca,'fontsize',16)
    xticks(inds1)
    xlim([inds1(1)-1,inds1(end)+1])
    ylim([.1,.9]);
    yticks(.1:.1:.9);
    
    subplot(1,3,2);hold on;
    errorbar(inds1,nanmean(ddanger(ii,i0:end),1),nanstd(ddanger(ii,i0:end),[],1)./sqrt(sum(1-isnan(ddanger(ii,i0:end)))),'color',cW,'marker','none','linestyle','none','LineWidth',0.75)
    errorbar(inds2,nanmean(ddanger(jj,i0:end),1),nanstd(ddanger(jj,i0:end),[],1)./sqrt(sum(1-isnan(ddanger(jj,i0:end)))),'color',cM,'marker','none','linestyle','none','linewidth',0.75)
    errorbar(inds3,nanmean(ddanger(kk,i0:end),1),nanstd(ddanger(kk,i0:end),[],1)./sqrt(sum(1-isnan(ddanger(kk,i0:end)))),'color',cS,'marker','none','linestyle','none','linewidth',0.75)
    plot(inds1,nanmean(ddanger(ii,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cW,'color',cW,'LineWidth',0.75)
    plot(inds2,nanmean(ddanger(jj,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cM,'color',cM,'LineWidth',0.75)
    plot(inds3,nanmean(ddanger(kk,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cS,'color',cS,'LineWidth',0.75)
    ylabel('distance to safety')
    xlabel('trial')
    set(gca,'fontsize',16)
    xticks(inds1)
    xlim([inds1(1)-1,inds1(end)+1])
    ylim([0,16]);
    yticks([0,4,8,12,16]);
    yticklabels({num2str(0),num2str(4*90./24),num2str(8*90./24),num2str(12*90./24),num2str(16*90./24)});
    
    subplot(1,3,3);hold on;
    errorbar(inds1,nanmean(PI(ii,i0:end),1),nanstd(PI(ii,i0:end),[],1)./sqrt(sum(1-isnan(PI(ii,i0:end)))),'color',cW,'marker','none','linestyle','none','LineWidth',0.75)
    errorbar(inds2,nanmean(PI(jj,i0:end),1),nanstd(PI(jj,i0:end),[],1)./sqrt(sum(1-isnan(PI(jj,i0:end)))),'color',cM,'marker','none','linestyle','none','linewidth',0.75)
    errorbar(inds3,nanmean(PI(kk,i0:end),1),nanstd(PI(kk,i0:end),[],1)./sqrt(sum(1-isnan(PI(kk,i0:end)))),'color',cS,'marker','none','linestyle','none','linewidth',0.75)
    plot(inds1,nanmean(PI(ii,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cW,'color',cW,'LineWidth',0.75)
    plot(inds2,nanmean(PI(jj,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cM,'color',cM,'LineWidth',0.75)
    plot(inds3,nanmean(PI(kk,i0:end),1),'o-','markersize',12,'markeredgecolor','none','markerfacecolor',cS,'color',cS,'LineWidth',0.75)
    ylabel('PI score')
    xlabel('trial')
    set(gca,'fontsize',16)
    xticks(inds1)
    xlim([inds1(1)-1,inds1(end)+1])
    legend({'','','','weak initial goal','medium initial goal', 'strong initial goal'});
end

if ismember(4,plotNum)

    pvals = [];

    rvec = [linspace(-1,1,24),linspace(1,-1,24),linspace(-1,1,24),linspace(1,-1,24)];
    for i=1:numel(sL.fly)
        for j=1:numel(sL.fly(1).trial)
            x  = rsL.fly(i).trial(j).x;
            xw = wrap(x,rsL.fly(i).punishedBar);
            
            inds = sL.fly(i).trial(j).PIinds; 
            xw = xw(inds);
            
            ddangerL(i,j)   = min(abs(sL.fly(i).trial(j).shiftedFull.prefLoc-[23,24,71,72]));
            goalLocL(i,j) = sL.fly(i).trial(j).shiftedFull.prefLoc*2*pi./96;        
            x  = quickWrap(xw + sL.fly(i).trial(j).shifted.shift,nx);
            rL(i,j) = (numel(find(x>=12 & x<36))-numel([find(x<12);find(x>=36)]))./numel(x);
        end
    end

    for i=1:numel(sKC.fly)
        for j=1:numel(sKC.fly(1).trial)
            x  = rsKC.fly(i).trial(j).x;
            xw = wrap(x,rsKC.fly(i).punishedBar);
            
            inds = sKC.fly(i).trial(j).PIinds; 
            xw = xw(inds);
            
            ddangerKC(i,j) = min(abs(sKC.fly(i).trial(j).shiftedFull.prefLoc-[23,24,71,72]));
            goalLocKC(i,j) = sKC.fly(i).trial(j).shiftedFull.prefLoc*2*pi./96;        
            x  = quickWrap(xw + sKC.fly(i).trial(j).shifted.shift,nx);
            rKC(i,j) = (numel(find(x>=12 & x<36))-numel([find(x<12);find(x>=36)]))./numel(x);
            PIKC(i,j)    = sKC.fly(i).trial(j).PI;
        end
    end



    %plot CDF of initial conditions
    figure; hold on;set(gcf,'Position',[200 200 400 800],'color','w')
    subplot(2,1,1);hold on;
    ii = find(~isnan(rL(:,2)));
    jj = find(~isnan(rKC(:,2)));
    plot(sort(rL( ii,2)),linspace(0,1,numel(ii)),'color','k','linewidth',2)
    plot(sort(rKC(jj,2)),linspace(0,1,numel(jj)),'color',[.6,.6,.6],'linewidth',2)
    xlim([-.2,1])
    ylim([0,1])
    xlabel('goal strength')
    ylabel('CDF')
    set(gca,'fontsize',16)
    legend({'WT','SS0009X'})

    subplot(2,1,2);hold on;
    ii = find(~isnan(goalLocL(:,2)));
    jj = find(~isnan(goalLocKC(:,2)));
    plot(sort(goalLocL( ii,2)),linspace(0,1,numel(ii)),'color','k','linewidth',2)
    plot(sort(goalLocKC(jj,2)),linspace(0,1,numel(jj)),'color',[.6,.6,.6],'linewidth',2)
    xlim([0,2*pi])
    ylim([0,1])
    xlabel('goal location')
    ylabel('CDF')
    set(gca,'fontsize',16)

end

end




