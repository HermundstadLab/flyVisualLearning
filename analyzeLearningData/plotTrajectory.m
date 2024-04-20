function plotTrajectory(rs,fly,trial)
% PLOTTRAJECTORY plots the unwrapped heading trajectories of an individual 
% fly, and superimposes this over safe/danger zones. 
%
% INPUTS:
%   rs: data structure containing the raw heading trajectories
%   fly: integer index of fly to plot
%   trial: integer index of trial to plot
%
% See also: PLOTSACCADES, PLOTBEHAVIORTRACES
%


nT = numel(rs.fly(1).trial(1).x);
t  = (1:nT)./1000;

x = unwrap(rs.fly(fly).trial(trial).x);
x1 = x;
if max(abs(diff(x)))>2
    ii = find(abs(diff(x))>2);
    if numel(ii)>2
        error('multiple wrapping errors')
        return
    end
end
x = x1;

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

cmap = parula(15);
cF = .5*cmap(1,:);
plot(t,x,'color',cF,'linewidth',3)
ylim([ymin,ymax])

