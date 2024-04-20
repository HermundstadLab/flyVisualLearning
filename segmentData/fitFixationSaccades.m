function fitFixationSaccades(sL,sNL,p)
% FITFIXATIONSACCADES fits parametric distributions to the duration of
% fixations and the duration and angular velocity of saccades, and saves
% the results
%
% INPUTS:
%   sL:  data structure containining segmented data for laser-trained flies
%   sNL: data structure containining segmented data for no-laser control
%       flies
%   p: parameter vector
%
% See also: GETDISTS

datasets = {sL,sNL};
nSets    = numel(datasets);
nFlies   = [numel(sL.fly),numel(sNL.fly)];
nTrials  = [numel(sL.fly(1).trial),numel(sNL.fly(1).trial)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           fit duration of fixations to individual trials/flies          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distribution = 'InverseGaussian';
threshold    = p.durThreshold;           

output  = getDists(datasets,'flies','fixations','dt',nTrials,nFlies);
muT      = nan(nSets,max(nFlies));
lambdaT  = nan(nSets,max(nFlies));
empMean = nan(nSets,max(nFlies));
empVar  = nan(nSets,max(nFlies));
for i=1:nSets
    for j=1:nFlies(i)
        M = output{i,j}';
        M(M<threshold) = [];
        d = fitdist(M,distribution);
        muT(     i,j) = d.mu;
        lambdaT( i,j) = d.lambda;
        empMean(i,j) = mean(M);
        empVar( i,j) = var( M);
    end
end
save('fits/fixDurFlies.mat','distribution','muT','lambdaT','empMean','empVar');


output  = getDists(datasets,'trials','fixations','dt',nTrials,nFlies);
muT      = nan(nSets,max(nTrials));
lambdaT  = nan(nSets,max(nTrials));
empMean = nan(nSets,max(nTrials));
empVar  = nan(nSets,max(nTrials));
for i=1:nSets
    for j=1:nTrials(i)
        M = output{i,j}';
        M(M<threshold) = [];
        d = fitdist(M,distribution);
        muT(    i,j) = d.mu;
        lambdaT(i,j) = d.lambda;
        empMean(i,j) = mean(M);
        empVar( i,j) = var( M);
    end
end
save('fits/fixDurTrials.mat','distribution','muT','lambdaT','empMean','empVar');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        fit angular velocity of saccades to individual trials/flies      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distribution = 'lognormal';
threshold    = p.velThreshold;

output  = getDists(datasets,'flies','saccades','wbaAvg',nTrials,nFlies);
muT      = nan(nSets,max(nFlies));
sigmaT   = nan(nSets,max(nFlies));
empMean = nan(nSets,max(nFlies));
empVar  = nan(nSets,max(nFlies));
for i=1:nSets
    for j=1:nFlies(i)
        M = abs(output{i,j})';
        M(M<threshold) = [];
        d = fitdist(M,distribution);
        muT(     i,j) = d.mu;
        sigmaT(  i,j) = d.sigma;
        empMean(i,j) = mean(M);
        empVar( i,j) = var( M);
    end
end
save('fits/sacVelFlies.mat','distribution','muT','sigmaT','empMean','empVar');


output  = getDists(datasets,'trials','saccades','wbaAvg',nTrials,nFlies);
muT      = nan(nSets,max(nTrials));
sigmaT   = nan(nSets,max(nTrials));
empMean = nan(nSets,max(nTrials));
empVar  = nan(nSets,max(nTrials));
for i=1:nSets
    for j=1:nTrials(i)
        M = abs(output{i,j})';
        M(M<threshold) = [];
        d = fitdist(M,distribution);
        muT(     i,j) = d.mu;
        sigmaT(  i,j) = d.sigma;
        empMean(i,j) = mean(M);
        empVar( i,j) = var( M);
    end
end
save('fits/sacVelTrials.mat','distribution','muT','sigmaT','empMean','empVar');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              fit duration of saccades conditioned on velocity           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------extract distribution of saccades-------------------%
distribution = 'InverseGaussian';
threshold    = p.velThreshold;

outputV  = getDists(datasets,'trials','saccades','wbaAvg',nTrials,nFlies);
outputT  = getDists(datasets,'trials','saccades','dt',nTrials,nFlies);

V = [];
T = [];
for i=1:nSets
    for j=1:nTrials(i)
        V = [V;abs(outputV{i,j})'];
        T = [T;outputT{i,j}'];
    end
end
ii = find(V<threshold);
V(ii) = [];
T(ii) = [];


%------------fit IG to condition distribution p( dt | v )--------------%

ddv = .05; 
vel = .25:ddv:2;
muT      = nan(size(vel));
lambdaT  = nan(size(vel));
empMean = nan(size(vel));
empVar  = nan(size(vel));

for i=1:numel(vel)
    jj=find(V>vel(i) & V<(vel(i)+ddv));
    d = fitdist(T(jj),distribution);
    
    muT(i)      = d.mu;
    lambdaT(i)  = d.lambda;
    empMean(i) = mean(T(jj));
    empVar(i)  = var( T(jj));
end

%--------fit sigmoid to relationship btw drift and saccade velocity-------%
m = log(muT.^3);
v = log(muT.^3./lambdaT);

b0 = mean(v-m);
a  = exp(-b0./2);
nu = a./muT;

x0 = [3,5,0,0];
sigmoidFit = @(x,xdata)(x(1)./(1+exp(-x(2).*(xdata-x(3)))) - x(4));
xlsq = lsqcurvefit(sigmoidFit,x0,vel,nu);


save('fits/sacDur.mat','distribution','vel','muT','lambdaT','empMean','empVar','sigmoidFit','xlsq','nu');



end

