function output = getDists(datasets,cut,type,field,nTrials,nFlies)
% GETDISTS reformats data to be fit to parametric distributions
%
% INPUTS:
%   datasets:  cell array containing datasets to be processed
%   cut: string that specifies whether to cut data by flies or trials. Can
%       take values of 'flies' or 'trials'
%   type: string that specifies which data type to fit to. Can take values
%       of 'fixations' or 'saccades'
%   field: string that specifies which data property to fit to. Can take
%       values of 'dt' or 'wbaAvg'
%   nTrials: array that specifies trial indicies
%   nFlies: array that specifies fly indicies
%
% See also: FITFIXATIONSACCADES 


%cut:   trials or flies
%type:  fixation or saccades
%field: duration, direction, vel

nSets = numel(datasets);
if strcmp(cut,'trials')
    nRuns1 = nTrials;
    nRuns2 = nFlies;
elseif strcmp(cut,'flies')
    nRuns1 = nFlies;
    nRuns2 = nTrials;
else
    error('unrecognized cut of data');
end
    

for i = 1:nSets
    rs = datasets{i};
    for j=1:nRuns1(i)
        
        M = [];
        for k=1:nRuns2(i)
            if strcmp(cut,'trials')
                if isfield(rs.fly(k).trial(j).(type),field)
                    M = [M,rs.fly(k).trial(j).(type).(field)];
                end
            else
                if isfield(rs.fly(j).trial(k).(type),field)
                    M = [M,rs.fly(j).trial(k).(type).(field)];
                end
            end
        end
        
        if strcmp(field,'dt') || strcmp(field,'dtDanger')
            M = M./1000;
        end
        
        output{i,j} = M;
        
    end
    
end

end