function plotEmpiricalPDFs(sL,sNL,p)
% PLOTEMPIRICALPDFS plots the empirical probability distributions of 
% saccade and fixation properties
%
% INPUTS:
%   sL, sNL:  data structures containining segmented data
%   p: structure containing parameter values
%

datasets = {sL,sNL};
nSets    = numel(datasets);
nFlies   = [44,40];
nTrials  = [9,9];

durThreshold = p.durThreshold;
velThreshold = p.velThreshold;

fixV  = getDists(datasets,'trials','fixations','wbaAvg',nTrials,nFlies);
fixT  = getDists(datasets,'trials','fixations','dt',nTrials,nFlies);
sacV  = getDists(datasets,'trials','saccades','wbaAvg',nTrials,nFlies);
sacT  = getDists(datasets,'trials','saccades','dt',nTrials,nFlies);

FV = [];
FT = [];
SV = [];
ST = [];
for i=1:nSets
    for j=1:nTrials(i)
        FV = [FV;fixV{i,j}'];
        FT = [FT;fixT{i,j}'];
        SV = [SV;sacV{i,j}'];
        ST = [ST;sacT{i,j}'];
    end
end
ii = find(FT<durThreshold);
jj = find(abs(SV)<velThreshold);

FV(ii) = [];
FT(ii) = [];
SV(jj) = [];
ST(jj) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          VELOCITY DISTRIBUTIONS                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             fixations: compute empirical PDF across all data            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;hold on;

db = .05;
bins = -1:db:1;
centers = (-1+db/2):db:1;
c = histcounts(FV,bins);
c = c./sum(c);

subplot(3,2,1);hold on;
stairs(centers(1:numel(c))-db/2,c,'color',[0,200,166]./255,'linewidth',1)

xlabel('angular velocity of fixations')
ylabel('PDF')
xlim([-2,2])
set(gca,'fontsize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              saccades: compute empirical PDF across all data            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,2,3);hold on;

db = .1;
bins = 0:db:5;
c = histcounts(SV(SV>0),bins);
c = c./sum(c);
stairs(bins(1:end-1),c,'color',[215,215,255]./255,'linewidth',1)

c = histcounts(-SV(SV<0),bins);
c = c./sum(c);
stairs(-fliplr(bins(2:end)),fliplr(c),'color',[215,215,255]./255,'linewidth',1)

fc = 100;
centers2 = db/2/fc:db/fc:5;
d = fitdist(SV(SV>0),'lognormal');
y = pdf('lognormal',centers2,d.mu,d.sigma);y = y./sum(y);
plot(centers2,y*fc,'color',[175,175,255]./255,'linewidth',1)

d = fitdist(-SV(SV<0),'lognormal');
y = pdf('lognormal',centers2,d.mu,d.sigma);y = y./sum(y);
plot(fliplr(-centers2),fliplr(y*fc),'color',[175,175,255]./255,'linewidth',1)

d = fitdist(abs(SV),'lognormal');
xlim([-2,2])
xlabel('angular velocity of saccades')
ylabel('PDF')
title(['\mu = ',num2str(round(d.mu,2)),'; \sigma = ',num2str(round(d.sigma,2))])
set(gca,'fontsize',16);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       saccades: plot directionality                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,5);hold on;

for k=1:2
    rs = datasets{k};
    for j=1:9
        frac=[];
        for i=1:nFlies(k)
            if isfield(rs.fly(i).trial(j).saccades,'dt')
                d = rs.fly(i).trial(j).saccades.dx';
                if numel(d)>0
                    frac = [frac;numel(d(d>0))./numel(d)];
                end

            end
        end
        D(k,j) = mean(frac);
    end
end


ma = {'o','s'};
D = 2*D-1;
for i=1:2
    scatter(D(i,:),normrnd(0,.01,[1,9]),200,'filled',ma{i},...
        'markerfacecolor','w','markeredgecolor',[175,175,255]./255,'linewidth',2)
end

ylim([-.5,.5])
xlim([-.2,.2])
yticks([]);yticklabels({});
xlabel('directional bias of saccades')
set(gca,'fontsize',16)
    
     
set(gcf,'Position',[200 200 800 1200],'color','w')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          DURATION DISTRIBUTIONS                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             fixations: compute empirical PDF across all data            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2);hold on;


db = .5;
bins = 0:db:100;
centers = db/2:db:100;
c = histcounts(FT,bins);
c = c./sum(c);

d = fitdist(FT,'inverseGaussian');
fc = 100;
centers2 = db/2/fc:db/fc:100;
y = pdf('inverseGaussian',centers2,d.mu,d.lambda);y = y./sum(y);

%plot empirical PDF and IG fit
stairs(centers(1:numel(c))-db/2,c,'color',[0,200,166]./255,'linewidth',1)
plot(centers2,y*fc,'color',[0,160,126]./255,'linewidth',2)


xlabel('fixation duration (s)')
ylabel('PDF')
xlim([0,20])
title(['median duration: ',num2str(round(median(FT),2)),'s'])
set(gca,'fontsize',16);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             saccades: compute empirical PDF across all data            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,2,4);hold on;

load('fits/sacDur.mat','vel','muT','lambdaT')
inds = [1,4,10];
v = abs(SV);
f = [.4,.7,1];

%plot empirical PDFs
db = .05;
bins = 0:db:5;
for i=1:numel(inds)
    kk=find( v>vel(inds(i)) & v<(vel(inds(i))+.05) );
    dt=ST(kk);

    c = histcounts(dt,bins);
    c = c./sum(c);
    stairs(bins(1:end-1),c,'color',f(i).*[175,175,255]./255,'linewidth',1)
end



fc = 100;

for m=1:3
    centers = db/2/fc:db/fc:5;
    y = pdf('inverseGaussian',centers,muT(inds(m)),lambdaT(inds(m)));y = y./sum(y);

    plot(centers,y*fc,'color',f(m).*[175,175,255]./255,'linewidth',2)

end
xlim([0,2])
xlabel('saccade duration (s)')
ylabel('PDF')
set(gca,'fontsize',16)
title(['median duration: ',num2str(round(median(ST),2)),'s'])

set(gcf,'Position',[200 200 1200 1000],'color','w')