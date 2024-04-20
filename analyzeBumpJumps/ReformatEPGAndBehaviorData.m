function procdata = ReformatEPGAndBehaviorData(res2PBeh, segmented2P)
% Take in background-subtracted fluorescence data from two-photon
% calcium imaging experiments from EPG neurons in the Ellipsoid Body (EB)
% of tethered flying flies. Compute DF/F and convert imaging and behavior
% data into a form that can be used to compute population vector
% averages, bump properties, bump jumps etc. The behavioral data is
% already segmented into saccades and fixations. The visual setting for all
% this data is 180-deg symmetric with four patterns that are 90-deg apart
% from each other, with matching pairs opposite each other.
%
% Note that this function moves arena x data so that the first "danger zone"
% is always in the middle. There is no danger zone here, because these
% experiments were performed without laser training, but the visual scene
% is the same one as used for laser-based operant learning in purely
% behavioral experiments that did involve such training.
%
% This uses the wrapFull function from AMH.
%
% Vivek Jayaraman

flydata = res2PBeh.fly;
procdata = flydata;

nFlies = length(flydata);
for f = 1:nFlies
    totNTrials = length(flydata(f).trial);
    setpointData = segmented2P.fly(f);

    for tr = 1:totNTrials
        tsteps = flydata(f).trial(tr).time; % imaging timestamps
        xpos_this = wrapFull(flydata(f).trial(tr).x, flydata(f).punishedBar); % Align to danger zone (punished bar)

        % The fly sometimes stops flying, or its wingbeat frequency drops below
        % acceptable levels. These timestamps are not considered for behavioral
        % analysis. Exclude for the rest of the analyses as well.
        PIinds = setpointData.trial(tr).PIinds;

        % Behavioral data is collected even after the imaging trial is done.
        % Remove those extra indices.
        PIinds = PIinds(PIinds <= length(xpos_this));
        tsteps = tsteps(PIinds);

        % Mark times when the fly is in safe and danger zones
        dangerVec_red = flydata(f).trial(tr).danger(PIinds);
        tDanger = length(find(dangerVec_red));
        tSafe = length(dangerVec_red) - tDanger;

        % Now compute a PI score that is based purely on this subset of data.
        % Note that this may differ from the PI score computed purely from
        % (longer) purely behavioral data
        PIscore = (tSafe - tDanger)/(tSafe+tDanger);

        % Fill analysis input data structure with appropriate time series data.
        procdata(f).trial(tr).x = xpos_this(PIinds);
        procdata(f).trial(tr).wba = flydata(f).trial(tr).wba(PIinds);
        procdata(f).trial(tr).wbf = flydata(f).trial(tr).wbf(PIinds);
        procdata(f).trial(tr).danger = dangerVec_red;
        procdata(f).trial(tr).puff = flydata(f).trial(tr).puff(PIinds);
        procdata(f).trial(tr).PIinds = PIinds;
        procdata(f).trial(tr).PI = PIscore; % Based purely on subset of indices here
        procdata(f).trial(tr).time = tsteps; % time
        lTrial = length(PIinds);

        % Copy over all saccade info as long as the saccades begin
        % within the bounds of the imaging part of the trial.
        procdata(f).trial(tr).saccades = [];
        if ~isempty(setpointData.trial(tr).saccades)
            saccadeIdxCondition = (setpointData.trial(tr).saccades.tstart < lTrial);
            if ~isempty(find(saccadeIdxCondition, 1))
                procdata(f).trial(tr).saccades = IndexedStructCopy(setpointData.trial(tr).saccades, saccadeIdxCondition);
            end
        end

        % Copy over all fixation info as long as the fixations begin
        % within the bounds of the imaging part of the trial.
        procdata(f).trial(tr).fixations = [];
        if ~isempty(setpointData.trial(tr).fixations)
            fixationIdxCondition = (setpointData.trial(tr).fixations.tstart < lTrial);
            if ~isempty(find(fixationIdxCondition, 1))
                procdata(f).trial(tr).fixations = IndexedStructCopy(setpointData.trial(tr).fixations, fixationIdxCondition);
            end
        end

        % Copy over the inferred preferred heading (based on
        % behavioral analysis, which may include more data than for imaging)
        procdata(f).trial(tr).setpoint = setpointData.trial(tr).shiftedFull.prefLoc;

        bgsubF = flydata(f).trial(tr).bgsubF;
        % And now, finally, tackle the imaging data.
        % F0: median of the lowest 10% of fluorescence across the trial.
        DF_F0 = median(mink(bgsubF, floor(0.1*size(bgsubF, 1)), 1), 1);
        DF_F0_mat = repmat(DF_F0, size(bgsubF, 1), 1);
        DFF_this = (bgsubF - DF_F0_mat)./DF_F0_mat;
        DFF_this = DFF_this(PIinds, :);
        DFF_this = DFF_this';

        procdata(f).trial(tr).dff = DFF_this(32:-1:1, :); % reordered
    end
end

end % of function

% Small function to copy over fields of a struct
function T = IndexedStructCopy(S, Condition, FieldList)
if nargin == 2
    FieldList = fieldnames(S);
end
for iField = 1:numel(FieldList)
    Field    = FieldList{iField};
    T.(Field) = S.(Field)(Condition);
end
end % of function