function [C,M,nOffsets,flyindex] = plotOffsetSummary(tSeriesProcAll,trials)
% PLOTOFFSETSUMMARY plots the number and angular separation of offsets
% between the bump in the EB and the visual scene. Summarized across flies
% and trials.
%
% INPUTS:
%   tSeriesProcAll: data structure containing fluorescence data
%   trials: integer indices of trials to use
%
% See also: PLOTOFFSETS
% 

M = nan(numel(tSeriesProcAll)*numel(trials),10);
C = [];
m = 1;
for fly=1:numel(tSeriesProcAll)
    
    D = nan(9,10);
    for i=1:numel(trials)
       
        trial = trials(i);
            
        if numel(tSeriesProcAll(fly).jumpData(trial).PIinds)>0
            density = tSeriesProcAll(fly).jumpData(trial).bumpjumps.density;
            xmesh   = tSeriesProcAll(fly).jumpData(trial).bumpjumps.xmesh;
            locs    = tSeriesProcAll(fly).jumpData(trial).bumpjumps.centerLocs;
            offsets = xmesh(locs);
            [~,minlocs] = findpeaks2(1-density);
            minlocs = [1,minlocs,numel(xmesh)];

            d = [];
            for j=1:numel(locs)
                i1 = find(minlocs<=locs(j),1,'last');
                i2 = find(minlocs>locs(j),1,'first');
                d(j) = sum(density(minlocs(i1):minlocs(i2)));
                if i1==1
                    d(j) = d(j)+sum(density(minlocs(end-1):minlocs(end)));
                end
                if i2==numel(xmesh)
                    d(j) = d(j)+sum(density(minlocs(1):minlocs(2)));
                end
            end
            d = d./sum(d);
            d(d<1./numel(xmesh)) = [];


            [~,isort] = sort(d,'descend');
            off_sort  = offsets(isort);

            M(m,1:numel(d)) = d(isort);
            flyindex(m,1:3) = [fly,mod(fly,2),trials(i)];
            D(i,1:numel(d)) = d(isort);

            if numel(off_sort)>1
                C(m) = min(abs(off_sort(1)-off_sort(2)),2*pi-abs(off_sort(1)-off_sort(2)));
            else
                C(m) = nan;
            end
            m = m+1;
        end

    end
    dmax = cumsum(nanmean(D));
    ii = find(dmax>.9,1,'first');
    N(fly) = ii;
            
end

theta = linspace(0,pi,17);
dtheta = theta(2)-theta(1);

cJump = histcounts(C,theta);
b = histcounts(N,0.5:1:5.5);


mMax = find(nansum(M,2)>eps,1,'last');
nMax = find(nansum(M,1)>eps,1,'last');


figure;set(gcf,'Position',[200 200 1000 400],'color','w')
Mbin = M;
Mbin(~isnan(M)) = 1;
bins = .5:1:5.5;
cOffset = histcounts(nansum(Mbin(1:mMax,:),2),bins);
subplot(1,2,1);bar(cOffset);
set(gca,'fontsize',16)
xlabel('number of offsets')
ylabel('count')

subplot(1,2,2);bar(theta(1:end-1)+dtheta./2,cJump)
xlim([0-.3,pi+.3])
xticks(0:pi/2:pi)
xticklabels({'0','\pi/2','\pi'})
set(gca,'fontsize',16)
xlabel('size of jumps (deg)')
ylabel('count')

B = nansum(Mbin,2);
for i=1:numel(tSeriesProcAll)
    for j=1:numel(trials)
        ii = find(flyindex(:,1)==i & flyindex(:,3)==j);
        if numel(ii)==1
            nOffsets(i,j) = B(ii);
        elseif numel(ii)>1
            disp('error');
        end
    end
end
