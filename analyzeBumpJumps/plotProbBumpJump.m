function [cJ,cV] = plotProbBumpJump(tSeriesProcAll,nOffsets,p)
% PLOTPROBBUMPJUMP computes and plots the conditional probability that the 
% compass bump in the EB jumps to a different orientation relative to the
% visual scene. This is measured relative to the inferred goal orientation 
% in the EB (i.e., the orientation at which the bump spends the most time 
% when the fly is orientating toward its preferred arena heading). The 
% conditional probability is computed for each flies, and then combined 
% across flies. The resulting histogram is then fit with a cosine.
%
% INPUTS:
%   tSeriesProcAll: data structure containing fluorescence data
%   nOffsets: array containing the number of unique offsets between the 
%       bump in the EB and the orientation of the visual scene, computed
%       across flies and trials
%   p: structure containing experimental parameters
%
% OUTPUTS:
%   cJ: number of bump jumps at different orientations in the EB, measured
%       relative to the inferred goal heading
%   cV: number of bump visits to different orientations in the EB, measured
%       relative to the inferred goal heading
%
% See also: PLOTJUMPS
%

cJ = [];cV = [];
for i=1:numel(tSeriesProcAll)
    [cVisit,cJump] = plotJumps(tSeriesProcAll,i,1:p.nTrials,nOffsets,p);
    close();
    cJ = [cJ;cJump];
    cV = [cV;cVisit];
end

theta = linspace(0,pi,17);
dtheta = theta(2)-theta(1);

%require 10 or more bump jumps
ii = find(sum(cJ,2)>9);
y = nanmean(cJ(ii,:)./cV(ii,:));
x = theta(1:end-1)+dtheta./2;
errfnc = @(c) sum((c(1).*cos(x)+c(2) - y).^2);
copt = fminunc(errfnc,[.2,.1]);


figure; hold on;hold on;set(gcf,'Position',[200 200 800 400],'color','w')
bar(x,nanmean(cJ(ii,:)./cV(ii,:)));plot(x,copt(1).*cos(x)+copt(2));
xlim([0-.3,pi+.3])
xticks(0:pi/2:pi)
xticklabels({'0','\pi/2','\pi'})
set(gca,'fontsize',16)
xlabel('distance to inferred goal heading')
ylabel('conditional probability of jump')