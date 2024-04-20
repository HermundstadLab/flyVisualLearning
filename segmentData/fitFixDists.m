function fitFixDists(sL,sNL)
% FITFIXDISTS checks the best-fitting parametric distribution for the
% duration of fixations, and saves the results
%
% INPUTS:
%   sL:  data structure containining segmented data for laser-trained flies
%   sNL: data structure containining segmented data for no-laser control
%       flies

%---------------remove fixations below a threshold------------------------%

%inverse gaussian is the best fitting distribution for a majority of
%conditions when the threshold is between 100 and 300. Above 300, the
%best-fitting distribution is a generalized pareto (commonly used to fit
%the tails of distributions). Below 100, the best-fitting distribution is
%lognormal

thresh = 200;

%---------------------- across individual trials -------------------------%

[distName,~,rankIG,rankIGBS] = fitParamDists_indTrials(sL,sNL,thresh);
save('fits/topRankingDists_indTrials.mat','distName','rankIG','rankIGBS')


%---------------------- across individual flies --------------------------%
[distName,~,rankIG,rankIGBS] = fitParamDists_indFlies(sL,sNL,thresh);
save('fits/topRankingDists_indFlies.mat','distName','rankIG','rankIGBS')

end


function [distName,params,rankIG,rankIGBS,BIC,empMean,empVar] = fitParamDists_indTrials(sL,sNL,thresh,dist)

if nargin<4
    dist = 'inverse gaussian';
end

s = {sL,sNL};
nFlies = [numel(sL.fly),numel(sNL.fly)];
nTrials = numel(sL.fly(1).trial);

params = zeros(2,nTrials,2);
rankIG = zeros(2,nTrials);
rankIGBS = zeros(2,nTrials);
BIC  = zeros(2,nTrials);
empMean  = zeros(2,nTrials);
empVar   = zeros(2,nTrials);
distName = cell(2,nTrials);

for i = 1:2
    rs = s{i};
    for j=1:9
        dt = [];
        for k=1:nFlies(i)
            if isfield(rs.fly(k).trial(j).fixations,'dt')
                dt = [dt,rs.fly(k).trial(j).fixations.dt];
            end
        end
        dt(dt<thresh)=[];
        dt = dt./1000;
        
        [d,~] = allfitdist(dt');
        
        for m=1:16
            if strcmp(d(m).DistName,dist)==1
                break
            end
        end
        
        for mm=1:16
            if strcmp(d(mm).DistName,'birnbaumsaunders')==1 || strcmp(d(mm).DistName,dist)==1
                break
            end
        end
        
        
        for p=1:numel(d(m).Params)
            params(i,j,p) = d(m).Params(p);
        end
        rankIG(i,j) = m;
        rankIGBS(i,j) = mm;
        distName{i,j} = d(1).DistName;
        BIC(i,j) = d(m).BIC;
        empMean(i,j) = mean(dt);
        empVar(i,j) = var(dt);
    end
    
end

end


function [distName,params,rankIG,rankIGBS,BIC,empMean,empVar] = fitParamDists_indFlies(sL,sNL,thresh,dist)

if nargin<5
    dist = 'inverse gaussian';
end

s = {sL,sNL};
nFlies = [numel(sL.fly),numel(sNL.fly)];
nFliesMax = max(nFlies);

params = zeros(2,nFliesMax,2);
rankIG = zeros(2,nFliesMax);
rankIGBS = zeros(2,nFliesMax);
BIC  = zeros(2,nFliesMax);
empMean  = zeros(2,nFliesMax);
empVar   = zeros(2,nFliesMax);
distName = cell(2,nFliesMax);

for i = 1:2
    rs = s{i};
    for j=1:nFlies(i)
        dt = [];
        for k=1:9
            if isfield(rs.fly(j).trial(k ).fixations,'dt')
                dt = [dt,rs.fly(j).trial(k).fixations.dt];
            end
        end
        dt(dt<thresh)=[];
        dt = dt./1000;
        
        [d,~] = allfitdist(dt');
        
        for m=1:16
            if strcmp(d(m).DistName,dist)==1
                break
            end
        end
        
        for mm=1:16
            if strcmp(d(mm).DistName,'birnbaumsaunders')==1 || strcmp(d(mm).DistName,dist)==1
                break
            end
        end
        
        for p=1:numel(d(m).Params)
            params(i,j,p) = d(m).Params(p);
        end
        rankIG(i,j) = m;
        rankIGBS(i,j) = mm;
        distName{i,j} = d(1).DistName;
        BIC(i,j) = d(m).BIC;
        empMean(i,j) = mean(dt);
        empVar(i,j) = var(dt);
        
    end
    
end

end
        