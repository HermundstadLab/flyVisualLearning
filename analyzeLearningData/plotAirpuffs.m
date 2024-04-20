function plotAirpuffs(rsL,sL,p)
% PLOTAIRPUFFS plots the distribution of airpuffs experienced by
% laser-trained and no-laser control flies
%
% INPUTS:
%   rsL: data structure containining raw data 
%   sL:  data structure containining segmented data
%   p: parameter vector
%

for i=1:numel(rsL.fly)
    for j=1:numel(rsL.fly(1).trial)
        x  = wrap(rsL.fly(i).trial(j).x,rsL.fly(i).punishedBar);
        ii = find(rsL.fly(i).trial(j).puff>0.5);
        dtPuff(i,j)   = numel(ii)./1000;
        xPuff{i,j}    = x(ii);

        PI(i,j)   = sL.fly(i).trial(j).PI;
        ninds(i,j) = numel(sL.fly(i).trial(j).PIinds);
    end
end
ii = find(sum(dtPuff,2)>0);
dtPuff = dtPuff(ii,:);
xPuff = xPuff(ii,:);

PInaive  = sum(ninds( ii,p.naiveTrials ).*PI(ii,p.naiveTrials ),2)./sum( ninds(ii,p.naiveTrials ),2);
PIprobe  = sum(ninds( ii,p.probeTrials2).*PI(ii,p.probeTrials2),2)./sum( ninds(ii,p.probeTrials2),2);
dPI = PIprobe - PInaive;


figure;set(gcf,'Position',[200 200 800 800],'color','w')
subplot(3,1,2); bar(sum(dtPuff,2));set(gca,'FontSize',16)
ylabel('total airpuff duration')

subplot(3,1,3); bar(dPI);set(gca,'FontSize',16)
ylim([-.4,.8])
yticks([-.4,0,.4,.8])
xlabel('flies')
ylabel('\Delta PI')

subplot(3,1,1);hold on;set(gca,'FontSize',16)
ylabel('location of puffs (deg)')
xlim([0,numel(ii)+1])
ylim([0,48])
yticks([0,12,24,36,48])
yticklabels({'-90','-45','0','45','90'})

cmap = parula(numel(ii));
nbins = 48;
bins = 0:nbins;
xpos = .5:1:(nbins-.5);
for i=1:numel(ii)
    xall = [];
    tall = [];
    trials = find(dtPuff(i,:)>0);
    for j=1:numel(trials)
        xall = [xall;xPuff{i,trials(j)}];
        tall = [tall;ones(numel(xPuff{i,trials(j)}),1)];
    end
    [~,~,bininds] = histcounts(xall,bins);
    dt = accumarray(bininds,tall,[nbins,1],@sum,0);

    xx = repmat(i,[nbins,1]);
    yy = xpos';
    zz = 200*dt/1000;
    jj = find(zz>0);
    scatter(xx(jj),yy(jj),zz(jj),cmap(i,:),'filled')
end


end