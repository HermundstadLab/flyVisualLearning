function plotSaccades(rs,s,fly,trial,pplot,plotSlopes)
% PLOTBEHAVIORTRACES plots the heading trajectories of an individual fly,
% segmented into fixations and saccades
%
% INPUTS:
%   rs: data structure containing the raw heading trajectories
%   s:  data structure containining segmented data
%   fly: integer index of fly to plot
%   trial: integer index of trial to plot
%   pplot: structure containing plotting parameters
%   plotSlopes: string that specifies whether to display slopes of
%       best-fitting lines fit to individual segments. Default is false.
%
% See also: PLOTTRAJECTORY, PLOTBEHAVIORTRACES
%

if nargin<6
    plotSlopes = false;
end


flyS = fly;
rsS  = s.fly(flyS).trial(trial).saccades;
rsF  = s.fly(flyS).trial(trial).fixations;
rsD  = s.fly(flyS).trial(trial).drift;
rsO  = s.fly(flyS).trial(trial).other;
nT = numel(rs.fly(flyS).trial(trial).x);
t  = (1:nT)./1000;

x = unwrap(rs.fly(fly).trial(trial).x);

dx = 48;
if strcmp(rs.fly(fly).punishedBar,'UpBar')
    xstart = 0;
else
    xstart = 24;
end
barStarts = [fliplr(xstart-dx:-dx:min(x)-dx),xstart:dx:max(x)+dx]; 
blankStarts = barStarts+24;
starts = sort([barStarts,blankStarts]);
ymin = starts(find(starts<min(x),1,'last'));
ymax = starts(find(starts>max(x),1,'first'))-1;


for i=1:numel(barStarts)
    fill([0,nT./1000,nT./1000,0],[barStarts(i),barStarts(i),barStarts(i)+23,barStarts(i)+23],[.9,.9,.9])
end

cS = [175,175,255]./255;
cF = pplot.cF;
cD = [.8,.8,.8];
cO = cD;
dt = 1;



iAll = [];
if isfield(rsO,'tstart')
    for i=1:numel(rsO.tstart)
        plot(t(rsO.tstart(i):dt:rsO.tend(i)),x(rsO.tstart(i):dt:rsO.tend(i)),'.','color',cO,'linewidth',6)
        iAll = [iAll,rsO.tstart(i):dt:rsO.tend(i)]; 
    end
end
if isfield(rsS,'tstart')
    for i=1:numel(rsS.tstart)
        plot(t(rsS.tstart(i):dt:rsS.tend(i)),x(rsS.tstart(i):dt:rsS.tend(i)),'.','color',cS,'linewidth',6)
        iAll = [iAll,rsS.tstart(i):dt:rsS.tend(i)]; 
    end
end
if isfield(rsF,'tstart')
    for i=1:numel(rsF.tstart)
        tfit = rsF.tstart(i):dt:rsF.tend(i);
        plot(t(tfit),x(tfit),'.','color',cF,'linewidth',6) 
        iAll = [iAll,rsF.tstart(i):dt:rsF.tend(i)]; 
        if plotSlopes
            yfit = rsF.linearSlope(i).*tfit+rsF.linearOffset(i);
            plot(t(tfit),yfit,'-k','linewidth',2)
        end
    end
end
if isfield(rsD,'tstart')
    for i=1:numel(rsD.tstart)
        tfit = rsD.tstart(i):dt:rsD.tend(i);
        plot(t(tfit),x(tfit),'.','color',cD,'linewidth',6) 
        iAll = [iAll,rsD.tstart(i):dt:rsD.tend(i)]; 
        if plotSlopes
            yfit = rsD.linearSlope(i).*tfit+rsD.linearOffset(i);
            plot(t(tfit),yfit,'-','color',[.6,.6,.6],'linewidth',2)
        end
    end
end

iAll = sort(iAll);
iRest = 1:dt:nT;
iRest(iAll) = [];
if numel(iRest)>0
    plot(t(iRest),x(iRest),'.r','linewidth',6)
end
ylim([ymin,ymax])
axis off
