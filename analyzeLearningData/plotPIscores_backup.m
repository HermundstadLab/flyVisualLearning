function [PIL,PINL,pvalsL,pvalsNL,pvals,nvalsL,nvalsNL] = plotPIscores(rsL,rsNL,sL,sNL,type,pattern,p,excludePuffs)
% PLOTPISCORES plots the evolution of PI scores across trials
%
% INPUTS:
%   rsL, rsNL: data structures containing the raw heading trajectories
%   sL, sNL:  data structures containining segmented data
%   type: string that specifies data type. Can take values 'laser' or 'kir'
%   pattern: string that specifies whether to align to preference or
%       safety. Can take values  of 'finalPref', 'trialPref', 'UpBar', 
%       'DnBar', 'both' (i.e., aligned to safety)
%   p: structure containing parameter values
%   excludePuffs: logical that specifies whether or not to combine exclude
%       flies that experienced airpuffs. Default is false.
%

if nargin<8
    excludePuffs = false;
end

%use wrapped data to compute PI scores
nx = p.nx/2;    

for i=1:numel(sL.fly)
    for j=1:numel(sL.fly(1).trial)
        if any(rsL.fly(i).trial(j).puff>.5)
            puffedL(i,j) = 1;
        else
            puffedL(i,j) = 0;
        end
        if strcmp(pattern,'finalPref')==1 || strcmp(pattern,'trialPref')==1
            x  = rsL.fly(i).trial(j).x;
            xw = wrap(x,rsL.fly(i).punishedBar);
            inds = sL.fly(i).trial(j).PIinds; 
            xw = xw(inds);
            
            if numel(inds)<p.PIindThreshold
                PIL(i,j) = nan;
            else
                if strcmp(pattern,'finalPref')==1
                    %shift wrt trial 9
                    x  = quickWrap(xw + sL.fly(i).trial(9).shifted.shift,nx);
                else 
                    %shift wrt current preference
                    x  = quickWrap(xw + sL.fly(i).trial(j).shifted.shift,nx);
                end
                PIL(i,j) = (numel(find(x>=nx/4 & x<3*nx/4))-numel([find(x<nx/4);find(x>=3*nx/4)]))./numel(x);
            end
        else
            PIL(i,j) = sL.fly(i).trial(j).PI; 
        end
        nindsL(i,j) = numel(sL.fly(i).trial(j).PIinds); 
    end
end

for i=1:numel(sNL.fly)
    for j=1:numel(sNL.fly(1).trial)
        if any(rsNL.fly(i).trial(j).puff>.5)
            puffedNL(i,j) = 1;
        else
            puffedNL(i,j) = 0;
        end
        if strcmp(pattern,'finalPref')==1 || strcmp(pattern,'trialPref')==1
            x  = rsNL.fly(i).trial(j).x;
            xw = wrap(x,rsNL.fly(i).punishedBar);
            inds = sNL.fly(i).trial(j).PIinds; 
            xw = xw(inds);

            if numel(inds)<p.PIindThreshold
                PINL(i,j) = nan;
            else            
                if strcmp(pattern,'finalPref')==1
                    %shift wrt trial 9
                    x  = quickWrap(xw + sNL.fly(i).trial(9).shifted.shift,nx);
                else 
                    %shift wrt current preference
                    x  = quickWrap(xw + sNL.fly(i).trial(j).shifted.shift,nx);
                end
                
                PINL(i,j) = (numel(find(x>=nx/4 & x<3*nx/4))-numel([find(x<nx/4);find(x>=3*nx/4)]))./numel(x);
            end
        else
            PINL(i,j) = sNL.fly(i).trial(j).PI; 
        end
        nindsNL(i,j) = numel(sNL.fly(i).trial(j).PIinds);       
    end
end


if strcmp(pattern,'UpBar')==1 || strcmp(pattern,'DnBar')==1
    for i=1:numel(rsL.fly)
        if strcmp(rsL.fly(i).punishedBar,pattern)==1
            barL(i) = 1;
        else
            barL(i) = 0;
        end
    end
    
    for i=1:numel(rsNL.fly)
        if strcmp(rsNL.fly(i).punishedBar,pattern)==1
            barNL(i) = 1;
        else
            barNL(i) = 0;
        end
    end
    
elseif strcmp(pattern,'both')==1 || strcmp(pattern,'finalPref')==1 || strcmp(pattern,'trialPref')==1
    barL( 1:numel( rsL.fly)) = 1;
    barNL(1:numel(rsNL.fly)) = 1;
else
    PILnaive = mean(PIL(:,p.naiveTrials),2);
    barL = zeros(numel(rsL.fly),1);
        
    PINLnaive = mean(PINL(:,p.naiveTrials),2);
    barNL = zeros(numel(rsNL.fly),1);
        
    if strcmp(pattern,'withPref')==1    
        barL(PILnaive>0) = 1;
        barNL(PINLnaive>0) = 1;
    elseif strcmp(pattern,'againstPref')==1
        barL(PILnaive<0) = 1;
        barNL(PINLnaive<0) = 1;
        
    else
        error('unrecognized pattern pref');
    end
end

%exclude flies that experienced airpuffs
if excludePuffs
    iremoveL  = find(sum(puffedL, 2)>0);
    iremoveNL = find(sum(puffedNL,2)>0);
    PIL( iremoveL, :) = [];
    PINL(iremoveNL,:) = [];
    nindsL( iremoveL, :) = [];
    nindsNL(iremoveNL,:) = [];
    barL( iremoveL, :) = [];
    barNL(iremoveNL,:) = [];
end

%remove

%plot results
figure;set(gcf,'Position',[200 200 1000 500],'color','w');hold on;
subplot(1,3,[1,2]);hold on;
bh=bar([nanmean(PIL,1)',nanmean(PINL,1)'],'facecolor','w');
bh(1).LineWidth = 2;
x0 = get(bh,'XData');
dx = get(bh,'XOffset');

errorbar(x0{1}+dx{1},nanmean(PIL,1),nanstd(PIL,[],1)./sqrt(sum(1-isnan(PIL))),'k','marker','none','linestyle','none','LineWidth',2)
errorbar(x0{2}+dx{2},nanmean(PINL,1),nanstd(PINL,[],1)./sqrt(sum(1-isnan(PINL))),'k','marker','none','linestyle','none')

xlabel('trial')
ylabel('PI score')
if strcmp(type,'laser')
    legend({'laser', 'no laser'},'Location','northwest')
elseif strcmp(type,'kir')
    legend({'parental controls','kir silenced'},'Location','northwest')
else
    error('unrecognized data type')
end
set(gca,'fontsize',16)
set(gcf,'Position',[200 200 1000 600])
ylim([-1,1])
for i=1:numel(sL.fly(1).trial)
    pL  = PIL( ~isnan(PIL( :,i)),i);
    pNL = PINL(~isnan(PINL(:,i)),i);
    pvalsL( i) = signrank(pL );
    pvalsNL(i) = signrank(pNL);
    pvals(  i) = ranksum(pL,pNL);
    nvalsL( i) = numel(pL );
    nvalsNL(i) = numel(pNL);    
end


PInaiveL  = sum(nindsL( :,p.naiveTrials).*PIL(:,p.naiveTrials),2)./sum( nindsL(:,p.naiveTrials),2);
PIprobe2L = sum(nindsL(:,p.probeTrials2).*PIL(:,p.probeTrials2),2)./sum(nindsL(:,p.probeTrials2),2);

PInaiveNL  = sum(nindsNL(:,p.naiveTrials).*PINL( :,p.naiveTrials),2)./sum( nindsNL(:,p.naiveTrials),2);
PIprobe2NL = sum(nindsNL(:,p.probeTrials2).*PINL(:,p.probeTrials2),2)./sum(nindsNL(:,p.probeTrials2),2);


PILall = [PInaiveL,PIprobe2L];
PILall(isnan(sum(PILall,2)),:) = [];
xL = repmat(4:5,[size(PILall,1 ),1]);
[~,ii] = sort(sum(PILall,2));
PILall = PILall(ii,:);

PINLall = [PInaiveNL,PIprobe2NL];
PINLall(isnan(sum(PINLall,2)),:) = [];
xNL = repmat(1:2,[size(PINLall,1 ),1]);
[~,ii] = sort(sum(PINLall,2));
PINLall = PINLall(ii,:);


p28L  = signrank(PILall(:,1),PILall(:,2));
p28NL = signrank(PINLall(:,1),PINLall(:,2));

p88 = ranksum(PILall(:,2),PINLall(:,2));
p22 = ranksum(PILall(:,1),PINLall(:,1));


subplot(1,3,3);hold on;

h = plot(xL',PILall','-o');
set(h, {'MarkerFaceColor'}, get(h,'Color'),'markeredgecolor','none'); 
title({['p_{NL} 2-8 = ',num2str(p28NL)];['p_L 2-8 = ',num2str(p28L)];['p_{L/NL} 2-2 = ',num2str(p22)];['p_{L/NL} 8-8 = ',num2str(p88)]});
set(gca,'fontsize',16)

h1 = plot(xNL',PINLall','-o');
set(h1, {'MarkerFaceColor'}, get(h1,'Color'),'markeredgecolor','none'); 
plot([0,6],[0,0],'--k')
xlim([0,6])
ylim([-1,1])
xticks([2,4])
xticklabels({'NL','L'})
set(gca,'fontsize',16)

end
