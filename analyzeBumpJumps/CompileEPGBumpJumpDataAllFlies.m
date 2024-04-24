function [tSeriesProcAll] = CompileEPGBumpJumpDataAllFlies(procData)
% CompileEPGBumpJumpDataAllFlies computes parameters associated with EPG
% calcium activity in the ellipsoid body and its relationship to visual
% stimulus and the fly's behavior. The main goal is to end up with a clear
% idea of where the EPG bump "jumps", i.e., when the offset between EPG
% activity and the fly's visual surroundings changes. This is usually
% because of visual scene symmetries at the local or global scale.
%
% INPUTS:
% procData(fly).trial(tr).x: position relative to arena reference;
% procData(fly).trial(tr).dff;
% procData(fly).trial(tr).wba;
% procData(fly).trial(tr).wbf;
% procData(fly).trial(tr).danger;
% procData(fly).trial(tr).puff;
%
% OUTPUTS:
% tSeriesAll: stores time series of behavior, stimulus, EPG activity,
%   etc.
%
% Vivek Jayaraman

%% Params for peak detection in calcium traces (using findpeaks)
threshPeakProminence = 1; % chosen based on std
threshPeakHeight = 0.5; % ignore peaks below this minimal threshold
threshPeakDist = 6; % ignore peaks separated by less than this
threshMaxPeakHeight = 1; % ignore if no peak is larger than this
wbaeps = 1; % threshold needed for wba
warning('off', 'signal:findpeaks:largeMinPeakHeight');

%% Basic setup of experiment
xpos0 = 0; dxpos = 2; xposN = 96; % arena params used to bin firing rate maps
binsHead = xpos0:dxpos:xposN;

% EB ROIs
nEBWedges = 32;
adjEBWedgeSteps = 2*pi/nEBWedges;
wedgeBearings = (((adjEBWedgeSteps/2)-pi):adjEBWedgeSteps:(pi-(adjEBWedgeSteps/2)))';
wedgeBearingsEdges = (-pi:adjEBWedgeSteps:pi)';

nFlies = length(procData);

for f = 1:nFlies
    nTrials = length(procData(f).trial);
    for tr = 1:nTrials
        if isempty(procData(f).trial(tr).PIinds)
            jumpData(tr).dffMat = [];
            jumpData(tr).arena = [];
            jumpData(tr).origPVA_EB = [];
            jumpData(tr).correctedPVA_EB = [];
            jumpData(tr).PVA_EB_strength = [];
            jumpData(tr).saccades = [];
            jumpData(tr).fixations = [];
            jumpData(tr).setpoint = [];
            jumpData(tr).PIinds = [];
            jumpData(tr).x = [];
            jumpData(tr).wbf = [];
            jumpData(tr).wba = [];
            jumpData(tr).time = [];
            jumpData(tr).trialType = [];
            jumpData(tr).trialStr = [];
            jumpData(tr).FRMaps = [];
            jumpData(tr).bumpjumps = [];
            continue;
        end
        data = procData(f).trial(tr);
        wba = data.wba; % WBA difference (for ~angular velocity)
        nTimeSteps = size(procData(f).trial(tr).dff, 2);
        tsTwoBumps = zeros(1, nTimeSteps);

        %% Sort out arena signal to world mapping
        % First the arena signal. Signal goes from 1 to 96 
        % (12 virtual panels).
        % We'll convert that to a circular world from -pi to pi.
        adjBearing = (data.x)*2*pi/96 - pi;
        adjBearing = adjBearing';

        % Now shift time so first time point is zero 
        data.time = data.time - data.time(1);

        xpos = data.x;
        
        % For firing rate maps/curves
        % Arena positions
        FRMapHeading = zeros(length(binsHead)-1, nEBWedges);
        MaxPeakLocMat = zeros(length(binsHead)-1, nEBWedges);
        FRMapHeadingCCW = zeros(length(binsHead)-1, nEBWedges);
        FRMapHeadingCW = zeros(length(binsHead)-1, nEBWedges);
        [NHead, HEdges, binH] = histcounts(xpos, binsHead);

        dffMat = data.dff(1:nEBWedges, :);

        %% Compute population vector average
        popVectAvgEBBearing = circ_mean(repmat(wedgeBearings,1,size(dffMat,2)), dffMat, 1);

        popVectAvgEBBearingAmp = circ_r(repmat(wedgeBearings,1,size(dffMat,2)), dffMat, [], 1);

        corrPopVectAvgEBBearing = [];

        actualPVAOffset = circ_dist(popVectAvgEBBearing, adjBearing);
        medianBearingOffset = circ_median(actualPVAOffset(1:100:end), 2);

        % Use the offset to correct predicted bearing by appropriate amount.
        corrPopVectAvgEBBearing = [corrPopVectAvgEBBearing circ_dist(popVectAvgEBBearing, ...
            repmat(medianBearingOffset, 1, size(popVectAvgEBBearing,2)))];

        % Initialize jumpdata struct 
        jumpData(tr).dffMat = dffMat;
        jumpData(tr).correctedPVA_EB = corrPopVectAvgEBBearing;
        jumpData(tr).arena = adjBearing;
        jumpData(tr).origPVA_EB = popVectAvgEBBearing;
        jumpData(tr).PVA_EB_strength = popVectAvgEBBearingAmp;
        jumpData(tr).saccades = data.saccades;
        jumpData(tr).fixations = data.fixations;
        jumpData(tr).setpoint = data.setpoint;
        jumpData(tr).PIinds = data.PIinds;
        jumpData(tr).x = data.x;
        jumpData(tr).wbf = data.wbf;
        jumpData(tr).wba = wba;
        jumpData(tr).time = data.time;
        jumpData(tr).trialType = data.trialType;
        jumpData(tr).trialStr = data.trialStr;

        [~, maxPeakWedges] = max(dffMat, [], 1, 'omitnan', 'linear');
        maxPeakWedgesMat = zeros(nEBWedges, length(xpos));
        maxPeakWedgesMat(maxPeakWedges) = 1;

        % For each EB ROI, compute mean firing rate across headings.
        for j = 1:nEBWedges
            meanDFF1DHeading_this = accumarray(binH, dffMat(j, :), [length(NHead) 1], @(x) mean(x));
            meanDFF1DHeading_this(isnan(meanDFF1DHeading_this)) = 0;
            FRMapHeading(:,j) = FRMapHeading(:,j) + meanDFF1DHeading_this;
        end

        % For each visual arena location, where's the preferred EB?
        for j = 1:nEBWedges
            peakLocInEB = accumarray(binH, maxPeakWedgesMat(j,:), [length(NHead) 1], @(x) sum(x));
            peakLocInEB(isnan(peakLocInEB)) = 0;
            MaxPeakLocMat(:,j) = MaxPeakLocMat(:,j) + peakLocInEB;
        end

        FRMaps.FRMapHeading = FRMapHeading;

        jumpData(tr).FRMaps = FRMaps;

        %% Jump detection
        % Detect multiple peaks in the bump profile, a classic sign of an
        % upcoming bump jump (although sometimes aborted)
        dffMatExt = [dffMat(:, :); dffMat(1:nEBWedges/2, :)];
        maxPeakVal = NaN*ones(size(popVectAvgEBBearing));
        maxPeak = NaN*ones(size(popVectAvgEBBearing));
        circDiffMaxPeaks = NaN*ones(size(popVectAvgEBBearing));

        for ts = 1:size(dffMat, 2)
            if (max(dffMat(:, ts)) > threshMaxPeakHeight)
                % for each time step, find locations of significant peaks.
                % The presence of two bumps may indicate an upcoming bump jump.
                [pks, locs] = findpeaks(smoothdata(dffMatExt(:, ts), 'sgolay', 6), 'MinPeakProminence', std(dffMatExt(:, ts))*threshPeakProminence, ...
                    'MinPeakDistance', threshPeakDist, 'MinPeakHeight', threshPeakHeight);
                if isempty(locs)
                    continue;
                end
                locs = locs(locs > threshPeakDist); % because anything within that distance could be a fake peak unless also captured by the repeating bit of the trace
                pks = pks(locs > threshPeakDist);

                locs = mod(locs, nEBWedges);
                [locs, uLocIdx] = unique(locs); pks = pks(uLocIdx);
                locs(locs == 0) = nEBWedges;

                [pks, idxsSort] = sort(pks, 'descend');
                locs = locs(idxsSort);
                if isempty(locs)
                    continue;
                end
                if (length(locs) > 1)  % If the largest peak is high enough and if there's more than one peak
                    if length(locs) > 2
                        locs = locs(1:2);
                        pks = pks(1:2);
                    end
                end

                maxPeak(ts) = wedgeBearings(locs(1));
                circDiffMaxPeaks(ts) = circ_dist(adjBearing(ts), maxPeak(ts));

                maxPeakVal(ts) = pks(1);
            end
        end

        % Find peaks in the bump offset distribution. Assign all bump offsets
        % to clusters associated with those peaks

        % Pad data to avoid any edge effects
        circDiffMaxPeaksChop = circDiffMaxPeaks(1:end);
        circDiffMaxPeaksExt = [circDiffMaxPeaksChop-2*pi,circDiffMaxPeaksChop,circDiffMaxPeaksChop+2*pi];
        circDiffMaxPeaksExt(circDiffMaxPeaksExt<-2*pi)=[];
        circDiffMaxPeaksExt(circDiffMaxPeaksExt>2*pi)=[];

        % Now look for peaks in bump offset distribution
        bins = linspace(-2*pi,2*pi,2*(32*96)+1);
        [~,~,ic] = histcounts(circDiffMaxPeaksExt,bins);
        ic(ic<1) = [];
        circDiffMaxPeaksDisc = bins(ic);

        % Get kernel density estimate
        if isMATLABReleaseOlderThan("R2023b")
            % Use Zdravko Botev's kde from matlabcentral
            [~, density, xmesh, ~] = kde_zb(circDiffMaxPeaksDisc,2^8);
        else
            [density, xmesh] = kde(circDiffMaxPeaksDisc,NumPoints=2^8);
        end

        % Find peaks in density distribution
        [~,locs]      = findpeaks(density);
       
        % Get troughs to use to make cuts for clusters that we'll use later.
        % Note that this cannot use findpeaks(), because we need to handle
        % cases when the 'peaks' are flat plateaus.
        locs_flip = find(islocalmax(1-density, 'MinProminence', 1e-4, ...
            'FlatSelection', 'center'));

        % Remove lower padding
        iremL = find(xmesh<-pi);
        locs(xmesh(locs)<-pi) = [];
        locs = locs-numel(iremL);
        locs_flip(xmesh(locs_flip)<-pi) = [];
        locs_flip = locs_flip-numel(iremL);
        density(xmesh<-pi)  = [];
        xmesh(xmesh<-pi)    = [];

        % Remove upper padding
        locs(xmesh(locs)>pi)=[];
        locs_flip(xmesh(locs_flip)>pi) = [];
        density(xmesh>pi) = [];
        xmesh(xmesh>pi)   = [];

        % Now classify all the offsets as belonging to whatever classes the
        % peaks of the distribution represent.
        % Data is circular, so need to take that into account when binning
        clusterCenters  = xmesh(locs);
        clusterCutCands = [-pi,xmesh(locs_flip),pi];
        if clusterCenters(1)>=clusterCutCands(1) && clusterCenters(1)<clusterCutCands(2)
            clusterCenters = [clusterCenters,clusterCenters(1)];
        else
            clusterCenters = [clusterCenters(end),clusterCenters];
        end
        [~,~,clusterIDs] = histcounts(circDiffMaxPeaks,clusterCutCands);

        % Grab nan points
        naninds  = find(clusterIDs<1);
        clusterIDs(naninds) = [];

        [clusterCenters,~,clusterIDs] = unique(clusterCenters(clusterIDs));
        centerLocs = locs;

        clusterIDs(clusterIDs<1) = [];
        [clusterCenters,~,clusterIDs] = unique(clusterCenters(clusterIDs));

        %% Save all kinds of parameters to tSeriesAll
        bumpjumps.xmesh = xmesh;
        bumpjumps.density = density;
        bumpjumps.centerLocs = centerLocs;
        bumpjumps.locsMaxSinglePeaks = maxPeak;
        bumpjumps.offsets = clusterCenters;
        jumpData(tr).bumpjumps = bumpjumps;
    end
    tSeriesProcAll(f).name = procData(f).name;
    tSeriesProcAll(f).jumpData = jumpData;
end

end % of function
