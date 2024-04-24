function plotOffsetsAcrossFlies(tSeriesProcAll,trials)
% PLOTOFFSETSACROSSFLIES plots the orientation and strength of different 
% offsets between the bump in the EB and the visual scene. Plotted for all 
% flies across multiple trials
%
% INPUTS:
%   tSeriesProcAll: data structure containing fluorescence data
%   trials:  integer indices of trials to use
%
% See also: PLOTOFFSETSUMMARY
%

figure;hold on;set(gcf,'Position',[200 200 1000 1000],'color','w')
xx = linspace(0,2*pi,100);

for fly=1:numel(tSeriesProcAll)
    m = 1;
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

            subplot(5,4,fly);hold on;
            r = 1;
            if m<2
                plot(r*cos(xx),r*sin(xx),'linewidth',.5,'color',[.8,.8,.8])
                plot(r/2*cos(xx),r/2*sin(xx),'linewidth',.5,'color',[.8,.8,.8])
                plot(r.*[-1,1],[0,0],'linewidth',.5,'color',[.8,.8,.8])
                plot([0,0],r.*[-1,1],'linewidth',.5,'color',[.8,.8,.8])
            end
            for j=1:numel(d)
                plot([0,d(j).*cos(offsets(j))],[0,d(j).*sin(offsets(j))],'-k','linewidth',2)
                scatter(d(j).*cos(offsets(j)),d(j).*sin(offsets(j)),400*d(j),d(j),'filled');
                if isMATLABReleaseOlderThan("R2022a")
                    caxis([0,1]);
                else
                    clim([0,1]);
                end
                colormap(flipud(gray));
            end
            pbaspect([1 1 1])
            axis off
            m = m+1;
        end
    end
end
