function [rsout] = extractSaccades(rs,rsout,fly,trial,tThresh,alignment)
% EXTRACTSACCADES partitions the heading trajectory of a fly into
% fixations and saccades
% 
% INPUTS:
%   rs: data structure containing the raw heading trajectories
%   rsout: data structure containining preprocessed data for segmentations 
%   fly: numeric index of fly to be segmented
%   trial: numeric index of trial to be segmented
%   arena: a string that specifies which flight arena is used. Can take 
%       values 'upright' or '2P'
%   alignment: a string that specifies how to align the behavioral data. 
%       Can take values 'safe', 'UpBar', 'DnBar', or 'none'. Default is
%       'safe'
%   plotOut: logical input to call plotting functions. Default is 'false'
%
% OUTPUT:
%   rsout: a data structure that contains segmentations into saccades and
%       fixations
%
% See also: EXTRACTSACCADESPREPROCESS_ALLFLIES

if nargin<6
    alignment = 'safe';
end



nT = numel(rs.fly(fly).trial(trial).x);
dx = 96;

if strcmp(alignment,'safe')==1
    xwrap     = wrap(rs.fly(fly).trial(trial).x,rs.fly(fly).punishedBar);
    xwrapFull = wrapFull(rs.fly(fly).trial(trial).x,rs.fly(fly).punishedBar);
elseif strcmp(alignment,'UpBar')==1
    xwrap     = wrap(rs.fly(fly).trial(trial).x,'DnBar');
    xwrapFull = wrapFull(rs.fly(fly).trial(trial).x,'DnBar');
elseif strcmp(alignment,'DnBar')==1
    xwrap     = wrap(rs.fly(fly).trial(trial).x,'UpBar');
    xwrapFull = wrapFull(rs.fly(fly).trial(trial).x,'UpBar');
elseif strcmp(alignment,'none')==1
    xwrap     = wrap(rs.fly(fly).trial(trial).x,'none');
    xwrapFull = wrapFull(rs.fly(fly).trial(trial).x,'none');
else
    error('unrecognized alignment')
end
x     = unwrap(rs.fly(fly).trial(trial).x);
wba   = rs.fly(fly).trial(trial).wba;
wbf   = rs.fly(fly).trial(trial).wbf;
wbaF  = bandpassmu(wba, 1000, 0.1, 10);
acc   = diff(wbaF);acc(nT) = acc(end);
laser = rs.fly(fly).trial(trial).laser;
tsafe = zeros(1,nT);
tsafe(xwrap>=12 & xwrap<36)=1; 
tlaser = zeros(1,nT);
tlaser(laser>.5)=1; 


%--------------- select segments of continous flight ---------------------%
tSaccade  = rsout.fly(fly).trial(trial).saccades.times0;
tFixation = rsout.fly(fly).trial(trial).fixations.times0;
sFixation = rsout.fly(fly).trial(trial).fixations.slopes0;
oFixation = rsout.fly(fly).trial(trial).fixations.offsets0;
tDrift    = rsout.fly(fly).trial(trial).drift.times0;
sDrift    = rsout.fly(fly).trial(trial).drift.slopes0;
oDrift    = rsout.fly(fly).trial(trial).drift.offsets0;
tOther    = rsout.fly(fly).trial(trial).other.times0;

%select segments of >=30s continuous flight
tFlightMin = 30000;

indsAll = [];
for count=1:size(tSaccade,1)
    indsAll = [indsAll,tSaccade(count,1):tSaccade(count,2)];
end

for count=1:size(tFixation,1)
    indsAll = [indsAll,tFixation(count,1):tFixation(count,2)];
end

for count=1:size(tDrift,1)
    indsAll = [indsAll,tDrift(count,1):tDrift(count,2)];
end

for count=1:size(tOther,1)
    indsAll = [indsAll,tOther(count,1):tOther(count,2)];
end

indsAll = sort(indsAll);

if numel(indsAll)>0
    indsKeep = [];
    
    segStarts = [indsAll(1),indsAll(find(diff(indsAll)>1)+1)];
    segStops =  [indsAll(diff(indsAll)>1),indsAll(end)];

    
    for i=1:numel(segStarts)
        nseg(i) = numel(segStarts(i):segStops(i));
    end
    
    %keep trial iff there is one segment that is at least 30sec long; 
    %within trial, keep all segments that exceed a duration tThresh
    if max(nseg)>=tFlightMin 
        for i=1:numel(segStarts)
            if numel(segStarts(i):segStops(i))>tThresh
                indsKeep = [indsKeep, segStarts(i):segStops(i)];
            end
        end
    end
else
    indsKeep = [];
end

%-------------------- rebuild saccade / fixation times -------------------%
tSaccadeKeep = [];
for count=1:size(tSaccade,1)
    stmp = tSaccade(count,1):tSaccade(count,2);
    if numel(intersect(stmp, indsKeep))==numel(stmp)
        tSaccadeKeep = [tSaccadeKeep;[tSaccade(count,1),tSaccade(count,2)]];
    end
end

tFixationKeep = [];
sFixationKeep = [];
oFixationKeep = [];
for count=1:size(tFixation,1)
    stmp = tFixation(count,1):tFixation(count,2);
    if numel(intersect(stmp, indsKeep))==numel(stmp)
        tFixationKeep = [tFixationKeep;[tFixation(count,1),tFixation(count,2)]];
        sFixationKeep = [sFixationKeep;sFixation(count,1)];
        oFixationKeep = [oFixationKeep;oFixation(count,1)];
    end
end

tDriftKeep = [];
sDriftKeep = [];
oDriftKeep = [];
for count=1:size(tDrift,1)
    stmp = tDrift(count,1):tDrift(count,2);
    if numel(intersect(stmp, indsKeep))==numel(stmp)
        tDriftKeep = [tDriftKeep;[tDrift(count,1),tDrift(count,2)]];
        sDriftKeep = [sDriftKeep;sDrift(count,1)];
        oDriftKeep = [oDriftKeep;oDrift(count,1)];
    end
end

tOtherKeep = [];
for count=1:size(tOther,1)
    stmp = tOther(count,1):tOther(count,2);
    if numel(intersect(stmp, indsKeep))==numel(stmp)
        tOtherKeep = [tOtherKeep;[tOther(count,1),tOther(count,2)]];
    end
end


%------------------------- store saccades --------------------------------%
if size(tSaccadeKeep,1)==0
    rsout.fly(fly).trial(trial).saccades.tstart  = [];
    rsout.fly(fly).trial(trial).saccades.tend    = [];
    rsout.fly(fly).trial(trial).saccades.dt      = [];
    rsout.fly(fly).trial(trial).saccades.xstart  = [];
    rsout.fly(fly).trial(trial).saccades.xend    = [];
    rsout.fly(fly).trial(trial).saccades.xstartF = [];
    rsout.fly(fly).trial(trial).saccades.xendF   = [];
    rsout.fly(fly).trial(trial).saccades.dx      = [];
    rsout.fly(fly).trial(trial).saccades.dxTot   = [];
    rsout.fly(fly).trial(trial).saccades.wbaMax  = [];
    rsout.fly(fly).trial(trial).saccades.wbaAvg  = [];
    rsout.fly(fly).trial(trial).saccades.wbfMax  = [];
    rsout.fly(fly).trial(trial).saccades.wbfAvg  = [];
    rsout.fly(fly).trial(trial).saccades.accMax  = [];
    rsout.fly(fly).trial(trial).saccades.accAvg  = [];
    rsout.fly(fly).trial(trial).saccades.dtSafe  = [];
    rsout.fly(fly).trial(trial).saccades.dtDanger        = [];
    rsout.fly(fly).trial(trial).saccades.startInDanger   = [];
    rsout.fly(fly).trial(trial).saccades.endInDanger     = [];
    rsout.fly(fly).trial(trial).saccades.startInLaser    = [];
    rsout.fly(fly).trial(trial).saccades.endInLaser      = [];
    rsout.fly(fly).trial(trial).saccades.timeSinceDanger = [];
    
else
    
    for count=1:size(tSaccadeKeep,1)
        iA = tSaccadeKeep(count,1);
        iB = tSaccadeKeep(count,2);

        rsout.fly(fly).trial(trial).saccades.tstart(count) = iA;
        rsout.fly(fly).trial(trial).saccades.tend(count)   = iB;
        rsout.fly(fly).trial(trial).saccades.dt(count)     = numel(iA:iB);
        rsout.fly(fly).trial(trial).saccades.xstart(count) = xwrap(iA);
        rsout.fly(fly).trial(trial).saccades.xend(count)   = xwrap(iB);
        
        rsout.fly(fly).trial(trial).saccades.xstartF(count) = xwrapFull(iA);
        rsout.fly(fly).trial(trial).saccades.xendF(count)   = xwrapFull(iB);
        
        rsout.fly(fly).trial(trial).saccades.dx(count)     = x(iB)-x(iA);  
        rsout.fly(fly).trial(trial).saccades.dxTot(count)  = sum(abs(diff(x(iA:iB))));  
        rsout.fly(fly).trial(trial).saccades.wbaMax(count) = max(wba(iA:iB)); 
        rsout.fly(fly).trial(trial).saccades.wbaAvg(count) = mean(wba(iA:iB));
        rsout.fly(fly).trial(trial).saccades.wbfMax(count) = max(wbf(iA:iB));
        rsout.fly(fly).trial(trial).saccades.wbfAvg(count) = mean(wbf(iA:iB));
        rsout.fly(fly).trial(trial).saccades.accMax(count) = max(acc(iA:iB));
        rsout.fly(fly).trial(trial).saccades.accAvg(count) = mean(acc(iA:iB));
        rsout.fly(fly).trial(trial).saccades.dtSafe(count) = numel(find(tsafe(iA:iB)>.5));
        rsout.fly(fly).trial(trial).saccades.dtDanger(count)      = numel(find(tsafe(iA:iB)<.5));
        rsout.fly(fly).trial(trial).saccades.startInDanger(count) = 1-tsafe(iA);
        rsout.fly(fly).trial(trial).saccades.endInDanger(count)   = 1-tsafe(iB);
        rsout.fly(fly).trial(trial).saccades.startInLaser(count)  = tlaser(iA);
        rsout.fly(fly).trial(trial).saccades.endInLaser(count)    = tlaser(iB);
        
        

    end
end
%-------------------------------------------------------------------------%


%-------------------- store fixation properties --------------------------%
if size(tFixationKeep,1)==0
    
    rsout.fly(fly).trial(trial).fixations.tstart   = [];
    rsout.fly(fly).trial(trial).fixations.tend     = [];
    rsout.fly(fly).trial(trial).fixations.dt       = [];
    rsout.fly(fly).trial(trial).fixations.dtDanger = [];
    rsout.fly(fly).trial(trial).fixations.dtLaser  = [];
    rsout.fly(fly).trial(trial).fixations.dxTot    = [];
    rsout.fly(fly).trial(trial).fixations.xavg     = [];
    rsout.fly(fly).trial(trial).fixations.xavgF    = [];
    rsout.fly(fly).trial(trial).fixations.xstd     = [];
    rsout.fly(fly).trial(trial).fixations.xrange   = [];
    rsout.fly(fly).trial(trial).fixations.wbaMax   = [];
    rsout.fly(fly).trial(trial).fixations.wbaAvg   = [];
    rsout.fly(fly).trial(trial).fixations.wbfMax   = [];
    rsout.fly(fly).trial(trial).fixations.wbfAvg   = [];
    rsout.fly(fly).trial(trial).fixations.wbfVar   = [];
    
else
    for count=1:size(tFixationKeep,1)
        iA = tFixationKeep(count,1);
        iB = tFixationKeep(count,2);

        xavg = median(x(iA:iB));
        xtest = [fliplr(xavg:-dx:min(x)-dx),xavg:dx:max(x)+dx];
        jj=find(xtest>=0 & xtest<dx,1,'first');

        if strcmp(alignment,'safe')==1
            xavg  =  wrap(xtest(jj),rs.fly(fly).punishedBar);
            xavgFull =  wrapFull(xtest(jj),rs.fly(fly).punishedBar);
        elseif strcmp(alignment,'UpBar')==1
            xavg  =  wrap(xtest(jj),'DnBar');
            xavgFull =  wrapFull(xtest(jj),'DnBar');
        elseif strcmp(alignment,'DnBar')==1
            xavg  =  wrap(xtest(jj),'UpBar');
            xavgFull =  wrapFull(xtest(jj),'UpBar');
        elseif strcmp(alignment,'none')==1
            xavg  =  wrap(xtest(jj),'none');
            xavgFull =  wrapFull(xtest(jj),'none');
        else
            error('unrecognized alignment')
        end


        rsout.fly(fly).trial(trial).fixations.tstart(count)   = iA;
        rsout.fly(fly).trial(trial).fixations.tend(count)     = iB;
        rsout.fly(fly).trial(trial).fixations.dt(count)       = numel(iA:iB);
        rsout.fly(fly).trial(trial).fixations.dtDanger(count) = numel(iA:iB)-numel(find(tsafe(iA:iB)));
        rsout.fly(fly).trial(trial).fixations.dtLaser(count)  = numel(find(tlaser(iA:iB)));
        rsout.fly(fly).trial(trial).fixations.dxTot(count)    = sum(abs(diff(x(iA:iB)))); 
        rsout.fly(fly).trial(trial).fixations.xavg(count)     = xavg;

        rsout.fly(fly).trial(trial).fixations.xavgF(count)  = xavgFull;
        rsout.fly(fly).trial(trial).fixations.xstd(count)   = std(x(iA:iB));
        rsout.fly(fly).trial(trial).fixations.xrange(count) = max(x(iA:iB))-min(x(iA:iB));
        rsout.fly(fly).trial(trial).fixations.wbaMax(count) = max(wba(iA:iB)); 
        rsout.fly(fly).trial(trial).fixations.wbaAvg(count) = mean(wba(iA:iB));
        rsout.fly(fly).trial(trial).fixations.wbfMax(count) = max(wbf(iA:iB)); 
        rsout.fly(fly).trial(trial).fixations.wbfAvg(count) = mean(wbf(iA:iB));
        rsout.fly(fly).trial(trial).fixations.wbfVar(count) = var(wbf(iA:iB));
        
        rsout.fly(fly).trial(trial).fixations.linearSlope( count) = sFixationKeep(count,:);
        rsout.fly(fly).trial(trial).fixations.linearOffset(count) = oFixationKeep(count,:);

    end
        
      
end
%-------------------------------------------------------------------------%

%-------------------- store drift properties --------------------------%
if size(tDriftKeep,1)==0
    
    rsout.fly(fly).trial(trial).drift.tstart   = [];
    rsout.fly(fly).trial(trial).drift.tend     = [];
    rsout.fly(fly).trial(trial).drift.dt       = [];
    rsout.fly(fly).trial(trial).drift.dtDanger = [];
    rsout.fly(fly).trial(trial).drift.dtLaser  = [];
    rsout.fly(fly).trial(trial).drift.dxTot    = [];
    rsout.fly(fly).trial(trial).drift.xavg     = [];
    rsout.fly(fly).trial(trial).drift.xavgF    = [];
    rsout.fly(fly).trial(trial).drift.xstd     = [];
    rsout.fly(fly).trial(trial).drift.xrange   = [];
    rsout.fly(fly).trial(trial).drift.wbaMax   = [];
    rsout.fly(fly).trial(trial).drift.wbaAvg   = [];
    rsout.fly(fly).trial(trial).drift.wbfMax   = [];
    rsout.fly(fly).trial(trial).drift.wbfAvg   = [];
    rsout.fly(fly).trial(trial).drift.wbfVar   = [];
    
else
    for count=1:size(tDriftKeep,1)
        iA = tDriftKeep(count,1);
        iB = tDriftKeep(count,2);

        xavg = median(x(iA:iB));
        xtest = [fliplr(xavg:-dx:min(x)-dx),xavg:dx:max(x)+dx];
        jj=find(xtest>=0 & xtest<dx,1,'first');

        if strcmp(alignment,'safe')==1
            xavg  =  wrap(xtest(jj),rs.fly(fly).punishedBar);
            xavgFull =  wrapFull(xtest(jj),rs.fly(fly).punishedBar);
        elseif strcmp(alignment,'UpBar')==1
            xavg  =  wrap(xtest(jj),'DnBar');
            xavgFull =  wrapFull(xtest(jj),'DnBar');
        elseif strcmp(alignment,'DnBar')==1
            xavg  =  wrap(xtest(jj),'UpBar');
            xavgFull =  wrapFull(xtest(jj),'UpBar');
        elseif strcmp(alignment,'none')==1
            xavg  =  wrap(xtest(jj),'none');
            xavgFull =  wrapFull(xtest(jj),'none');
        else
            error('unrecognized alignment')
        end


        rsout.fly(fly).trial(trial).drift.tstart(count)   = iA;
        rsout.fly(fly).trial(trial).drift.tend(count)     = iB;
        rsout.fly(fly).trial(trial).drift.dt(count)       = numel(iA:iB);
        rsout.fly(fly).trial(trial).drift.dtDanger(count) = numel(iA:iB)-numel(find(tsafe(iA:iB)));
        rsout.fly(fly).trial(trial).drift.dtLaser(count)  = numel(find(tlaser(iA:iB)));
        rsout.fly(fly).trial(trial).drift.dxTot(count)    = sum(abs(diff(x(iA:iB)))); 
        rsout.fly(fly).trial(trial).drift.xavg(count)     = xavg;

        rsout.fly(fly).trial(trial).drift.xavgF(count)  = xavgFull;
        rsout.fly(fly).trial(trial).drift.xstd(count)   = std(x(iA:iB));
        rsout.fly(fly).trial(trial).drift.xrange(count) = max(x(iA:iB))-min(x(iA:iB));
        rsout.fly(fly).trial(trial).drift.wbaMax(count) = max(wba(iA:iB)); 
        rsout.fly(fly).trial(trial).drift.wbaAvg(count) = mean(wba(iA:iB));
        rsout.fly(fly).trial(trial).drift.wbfMax(count) = max(wbf(iA:iB)); 
        rsout.fly(fly).trial(trial).drift.wbfAvg(count) = mean(wbf(iA:iB));
        rsout.fly(fly).trial(trial).drift.wbfVar(count) = var(wbf(iA:iB));
        
        rsout.fly(fly).trial(trial).drift.linearSlope( count) = sDriftKeep(count,:);
        rsout.fly(fly).trial(trial).drift.linearOffset(count) = oDriftKeep(count,:);
        
    end
        
      
end
%-------------------------------------------------------------------------%


%-------------------- store other properties --------------------------%
if size(tOtherKeep,1)==0
    
    rsout.fly(fly).trial(trial).other.tstart   = [];
    rsout.fly(fly).trial(trial).other.tend     = [];

else
    for count=1:size(tOtherKeep,1)
        iA = tOtherKeep(count,1);
        iB = tOtherKeep(count,2);

        rsout.fly(fly).trial(trial).other.tstart(count)   = iA;
        rsout.fly(fly).trial(trial).other.tend(count)     = iB;
    end
        
      
end
%-------------------------------------------------------------------------%



xwrap     = xwrap(indsKeep);
laser     = laser(indsKeep);
indDanger = find(xwrap<12 | xwrap>=36);
indSafe   = find(xwrap>=12 & xwrap<36);
PI        = (numel(indSafe)-numel(indDanger))./numel(xwrap);
nlaser    = numel(find(laser>.5));

rsout.fly(fly).trial(trial).PI     = PI;
rsout.fly(fly).trial(trial).PIinds = indsKeep;  
rsout.fly(fly).trial(trial).nlaser = nlaser;



end


function [turnsAccum,turnsRest] = cutTurns(turns,s)

turnsAccum = {};
turnsRest = {};
m = 1;
p = 1;

%break turns?
for i=1:numel(turns)
    ss=s(turns{i});
   
    inds=find(ss);
    if numel(inds)>2
        if mean(ss(inds))>0
            [~,istart] = max((diff([0;ss(inds)])));      
            iend = find(ss(inds(istart+1:end))<ss(inds(istart))/4,1,'first');
        else
            [~,istart] = max(-(diff([0;ss(inds)])));
            iend = find(ss(inds(istart+1:end))>ss(inds(istart))/4,1,'first');
        end

        if numel(iend)>0        
            if inds(istart)>1
                turnsAccum{m} = turns{i}(1):turns{i}(inds(istart)-1);
                turnsAccum{m+1} = turns{i}(inds(istart)):turns{i}(inds(istart+iend));
                if inds(iend)<numel(turns{i})
                    turnsRest{p} = turns{i}(inds(istart+iend)+1:end);
                    p=p+1;
                end
                ntot = numel(1:inds(istart)-1)+numel(inds(istart):inds(istart+iend))+numel(inds(istart+iend)+1:numel(turns{i}));
                if numel(turns{i})~=ntot
                    disp(['error: ',num2str([numel(turns{i}),ntot])])
                end
                m=m+2;
            else
                %disp(inds(istart))
                turnsAccum{m} = turns{i}(inds(istart)):turns{i}(inds(istart+iend));
                if inds(iend)<numel(turns{i})
                    turnsRest{p} = turns{i}(inds(istart+iend)+1:end);
                    p=p+1;
                end
                ntot = numel(inds(istart):inds(istart+iend))+numel(inds(istart+iend)+1:numel(turns{i}));
                if numel(turns{i})~=ntot
                    disp(['error: ',num2str([numel(turns{i}),ntot])])
                end
                m=m+1;
            end
        else
            turnsAccum{m} = turns{i};
            m=m+1;
        end
    else
        turnsAccum{m} = turns{i};
        m=m+1;
    end
end

end
function t = getTimes(c)
t = [];
for i=1:numel(c)
    t = [t,c{i}];
end
end

function [mc] = medcouple(x)
%
% 'medcouple' computes the medcouple measure, a robust measure of skewness
% for a skewed distribution. It takes into account cases where the
% observations are equal to the median of the series.
%
% Data in 'x' are organized so that columns are the time series and rows
% are the time intervals. All series contain the same number of
% observations.
%
% [mc] = medcouple(x) returns the following:
% mc    - vector with the medcouple measure of the data series
%
% Created by Francisco Augusto Alcaraz Garcia
%            alcaraz_garcia@yahoo.com
%
% References:
%
% 1) G. Brys; M. Hubert; A. Struyf (2004). A Robust Measure of Skewness.
% Journal of Computational and Graphical Statistics 13(4), 996-1017.

[n, c] = size(x);

[s_x,~] = sort(x);
x_med = nanmedian(x);
z = s_x - repmat(x_med,n,1);

mc = zeros(1,c);

for w = 1:c
    [ip, ~] = find(z(:,w)>=0); % These are the positions in z of z+
    [im, ~] = find(z(:,w)<=0); % These are the positions in z of z-

    p = size(ip,1);
    q = size(im,1);

    [mi, mj] = ind2sub([p,q],1:p*q); % Positions of all combinations of z+ and z- as elements in a pxq matrix

    zp = z(ip,w); % z+ repeated to account for all cells in the matrix
    zm = z(im,w); % z- repeated to account for all cells in the matrix

    h = (zp(mi)+zm(mj))./(zp(mi)-zm(mj)); % same size as mi, mj

    [ipz,~]= find(zp==0);   % row numbers of z+ = 0, i.e., x_{i} = median(x)
    [imz,~]= find(zm==0);   % row numbers of z- = 0, i.e., x_{i} = median(x)
    piz = ismember(mi,ipz); % positions in mi where z+=0
    pjz = ismember(mj,imz); % positions in mi where z-=0
    zmember = piz+pjz;      % same size as mi, mj
    pijz = find(zmember == 2);          % positions where z+ = z- = 0, i.e., x_{i} = x_{j} = median(x)
    [indi,indj] = ind2sub([p,q],pijz);  % pxq matrix position of the zero entries
    indi = indi - min(indi) + 1;        % row position of the zero entries as if they were in a separated matrix
    indj = indj - min(indj) + 1;        % column position of the zero entries as if they were in a separated matrix

    for i=1:size(pijz,2)
        if (indi(i) + indj(i) - 1) > size(find(z==0),1)
            h(pijz(i)) = 1;
        elseif (indi(i) + indj(i) - 1) < size(find(z==0),1)
            h(pijz(i)) = -1;
        else
            h(pijz(i)) = 0;
        end
    end
    mc(w) = median(h);
end

end


function turns = segTurns(fWBA)

nT = numel(fWBA);

changepts = find(diff(sign(fWBA))~=0);

starts = [1;changepts+1];
stops  = [changepts;nT];

turns = {};
for i=1:numel(starts)
    turns{i} = starts(i):stops(i);
end
end

