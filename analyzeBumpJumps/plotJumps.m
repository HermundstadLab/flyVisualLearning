function [cVisit,cJump] = plotJumps(tSeriesProcAll,fly,trials,nOffsets,p)
% PLOTJUMPS plots the orientation of the bump jumps in the EB for a single
% fly, accumulated across trials
%
% INPUTS:
%   tSeriesProcAll: data structure containing fluorescence data
%   fly:   integer index of fly to plot
%   trials:  integer indices of trials to use
%   nOffsets: array containing the number of unique offsets between the 
%       bump in the EB and the orientation of the visual scene, computed
%       across flies and trials
%   p: structure containing experimental parameters
%
% OUTPUTS:
%   cJump: histogram counts of number of bump jumps at different 
%       orientations in the EB, measured relative to the inferred goal 
%       heading
%   cVisit: histogram counts of number of bump visits to different 
%       orientations in the EB, measured    relative to the inferred goal 
%       heading
%  
% See also: PLOTPROBBUMPJUMP
%

jumpLocsEB = [];
jumpLocsArena = [];
dJumps = [];
sJumps = [];
hBall  = [];
hAall  = [];
hSPall = [];

theta = linspace(-pi,pi,33);
figure;hold on;
subplot(4,1,[1,2,3]);hold on;
for j=1:numel(trials)
    if numel(tSeriesProcAll(fly).jumpData(trials(j)).PIinds)>0 && nOffsets(fly,j)==2
        
        sp = tSeriesProcAll(fly).jumpData(trials(j)).setpoint*2*pi./p.nx;
        hBump  = tSeriesProcAll(fly).jumpData(trials(j)).bumpjumps.locsMaxSinglePeaks;
        hArena = tSeriesProcAll(fly).jumpData(trials(j)).arena;

        hB = hBump(1);
        hA = hArena(1);
        shifts = [-4*pi,-2*pi,0,2*pi,4*pi];
        tJumps = [];
        
        
        density = tSeriesProcAll(fly).jumpData(trials(j)).bumpjumps.density;
        xmesh   = tSeriesProcAll(fly).jumpData(trials(j)).bumpjumps.xmesh;
        locs    = tSeriesProcAll(fly).jumpData(trials(j)).bumpjumps.centerLocs;
        offsets = xmesh(locs);
        [~,minlocs] = findpeaks2(1-density);
        minlocs = [1,minlocs,numel(xmesh)];

        d = [];
        for jj=1:numel(locs)
            i1 = find(minlocs<=locs(jj),1,'last');
            i2 = find(minlocs>locs(jj),1,'first');
            d(jj) = sum(density(minlocs(i1):minlocs(i2)));
            if i1==1
                d(jj) = d(jj)+sum(density(minlocs(end-1):minlocs(end)));
            end
            if i2==numel(xmesh)
                d(jj) = d(jj)+sum(density(minlocs(1):minlocs(2)));
            end
        end
        d = d./sum(d);
        d(d<1./numel(xmesh)) = [];
        [~,isort] = sort(d,'descend');
        sp_shift = mod(sp-offsets(isort(1)),2*pi)-pi;

        hSP  = sp_shift;    
            
            
        for i=1:numel(hBump)-1
            x1 = hB(i);
            y1 = hA(i);
            x2 = repmat(hBump(i+1) +shifts,[1,numel(shifts)]);
            y2 = repmat((hArena(i+1)+shifts)',[1,numel(shifts)])';
            y2 = y2(:)';
            d = sqrt((x1-x2).^2 + (y1-y2).^2);
            [~,ii] = min(d);
            dmin = d(ii);
            dx   = abs(x2(ii)-x1);
            dy   = abs(y2(ii)-y1);
            hB   = [hB,x2(ii)];
            hA   = [hA,y2(ii)];    
            hSP  = [hSP,sp_shift];

            if dmin>3*pi/4 && dmin<5*pi/4
                tJumps = [tJumps,i];
                dJumps = [dJumps,dmin];
                
            end
        end



        for i=7*pi:-2*pi:pi
            hB(hB>i) = hB(hB>i)-2*pi;
        end
        for i=-7*pi:2*pi:-pi
            hB(hB<i) = hB(hB<i)+2*pi;
        end

        for i=7*pi:-2*pi:pi
            hA(hA>i) = hA(hA>i)-2*pi;
        end
        for i=-7*pi:2*pi:-pi
            hA(hA<i) = hA(hA<i)+2*pi;
        end

        for i=1:numel(tJumps)
            %only keep jumps if arena position didn't change instantaneously
            if abs(hA(tJumps(i))-hA(tJumps(i)+1))<1/p.nx
                jumpLocsEB    = [jumpLocsEB;   [hB(tJumps(i)),hB(tJumps(i)+1)]];
                jumpLocsArena = [jumpLocsArena;[hA(tJumps(i)),hA(tJumps(i)+1)]];
                sJumps = [sJumps;sp_shift];
            end
        end

        scatter(-hA,hB,50,[.8,.8,.8],'filled')
        hBall = [hBall,hB];
        hAall = [hAall,hA];
        hSPall = [hSPall,hSP];
    end
end


plot(-jumpLocsArena',jumpLocsEB','-k')
if numel(jumpLocsEB)>0
    scatter(-jumpLocsArena(:,1),jumpLocsEB(:,1),30,[.2,.2,.2],'filled')
    scatter(-jumpLocsArena(:,2),jumpLocsEB(:,2),80,[1 , 0, 0],'filled')
end
xlim([-pi-.3,pi+.3])
xticks(-pi:pi/2:pi)
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ylim([-pi-.3,pi+.3])
yticks(-pi:pi/2:pi)
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ylabel('compass heading')
xlabel('arena heading')
set(gca,'fontsize',16)


theta0 = linspace(-pi,pi,33);
dtheta0 = theta0(2)-theta0(1);

theta = linspace(0,pi,17);
dtheta = theta(2)-theta(1);
cDelta    = histcounts(dJumps,theta);

if numel(jumpLocsEB)>0
    
    distSPjump = min([abs(jumpLocsEB(:,1)-sJumps),2*pi-abs(jumpLocsEB(:,1)-sJumps)],[],2);
    distSPhead = min([abs(hBall'-hSPall'),2*pi-abs(hBall'-hSPall')],[],2);
    cJump  = histcounts(distSPjump,theta);
    cVisit = histcounts(distSPhead,theta);
    
    subplot(4,1,4);hold on;
    bar(theta(1:end-1)+dtheta./2,cJump./cVisit)
    set(gca,'fontsize',16)
    xlim([0-.3,pi+.3])
    xticks(0:pi/2:pi)
    xticklabels({'0','\pi/2','\pi'})
    xlabel('bump heading')
    

else
    cJump  = nan(1,16);
    cVisit = nan(1,16);
end