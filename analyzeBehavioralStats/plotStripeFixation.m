function [orientationData,strengthData] = plotStripeFixation(rsK90C,rsK90,rsK96C,rsK96,rsL,rsNL)
% PLOTSTRIPEFIXATION plots the circular mean of flies' heading trajectories
% during an initial stripe fixation trial.
%
% INPUTS:
%   rsK90C,rsK90,rsK96C,rsK96,rsL,rsNL: data structures containing raw 
%       heading trajectories
%
% OUTPUTS:
%   orientationData: cell array with orientations of preferred headings
%   strengthData:    cell array with strengths of preferred headings
%

rsAll = {rsK90C,rsK90,rsK96C,rsK96,rsL,rsNL};
titles = {'SS00090 Parent','SS00090 Kir','SS00096 Parent','SS00096 Kir','WT laser','WT no-laser'};
Tmax = 30000;
figure;set(gcf,'Position',[200 200 1400 600],'color','w');hold on;
for i=1:numel(rsAll)
    rs = rsAll{i};
    orientation = [];
    strength    = [];
    for j=1:numel(rs.fly)
        theta = rs.fly(j).stripeFixation.x*2*pi/96-pi;
        theta = theta(1:Tmax);
        orientation(j) = circ_mean(theta);
        strength(j)   = circ_r(theta);

        subplot(2,6,i);polarplot([0,orientation(j)],[0,strength(j)]);hold on
        rlim([0,1])
        title(titles{i});
        set(gca,'fontsize',16)

        subplot(2,6,i+6);plot(strength(j),1+normrnd(0,.1),'o','MarkerFaceColor','k','MarkerEdgeColor','none');hold on
    end
    orientationData{i} = orientation;
    strengthData{i} = strength;

    range = quantile(strength,[.25, .75]);
    subplot(2,6,i+6);
    plot(range,[1,1],'r','linewidth',2)
    plot(nanmedian(strength),1,'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',12)
    ylim([0,2])
    yticks([])
    xlabel('strength')
    set(gca,'fontsize',16)

end

