function tThresh = computeSegmentationThresh(sNL)
% COMPUTESEGMENTATIONTHRESH extracts the distribution of fixation durations
% across flies and trials within a dataset, and computed 75th quantile of
% the distribution
% 
% INPUTS:
%   sNL:  data structure containining segmented data from no-laser control
%       flies
%
% OUTPUT:
%   tThresh: 75th quantile of the distribution of fixation durations 
%


durations = [];
ntrials = 0;
for i=1:numel(sNL.fly)
    for j=1:numel(sNL.fly(i).trial)
        dtF = diff(sNL.fly(i).trial(j).fixations.times0,[],2)+1;
        durations = [durations; dtF];
        ntrials = ntrials+1;
    end
end
tThresh = quantile(durations,.75);

