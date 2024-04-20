function rsout = extractSaccadesPreprocess(rs,rsout,fly,trial,arena,alignment,plotOutput)
% EXTRACTSACCADESPREPROCESS partitions the heading trajectory of a fly into
% fixations and saccades
% 
% INPUTS:
%   rs: data structure containing the raw heading trajectories
%   rsout: data structure containining preprocessed data from other flies;
%       can be empty
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

if nargin<7
    plotOutput = false;
    if nargin<6
        alignment = 'safe';
    end
end

if strcmp(arena,'upright')
    wbfThresh = 0.01; 
elseif strcmp(arena,'2P')
    wbfThresh = 0.1;
else
    error('unrecognized arena input')
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
puff  = rs.fly(fly).trial(trial).puff;
tsafe = zeros(1,nT);
tsafe(xwrap>=12 & xwrap<36)=1; 
tlaser = zeros(1,nT);
tlaser(laser>.5)=1; 

%flag times before/after airpuff
pmin = 500;
pmax = 500;
indsPuff = find(puff>.5);
puffInds = [];
if numel(indsPuff)>0
    puffStarts = [indsPuff(1);indsPuff(find(diff(indsPuff)>1)+1)];
    puffStops  = [indsPuff(diff(indsPuff)>1);indsPuff(end)];
    for j=1:numel(puffStarts)
        if j>1
            psa = max([1,puffStops(j-1)+1,puffStarts(j)-pmin]);
        else
            psa = max([1,puffStarts(j)-pmin]);
        end
        if j<numel(puffStarts)
        	pst = min([120000,puffStarts(j+1)+1,puffStops(j)+pmax]);
        else
            pst = min([120000,puffStops(j)+pmax]);
        end
        puffInds = [puffInds,psa:pst];
    end
end

%compute changes in WBA that result in a change in position (used later to identify saccades)
s = -diff(x).*wbaF(1:end-1).*sign(diff(x));

%segment all turns
turns = segTurns(wbaF);                 
turns{end} = turns{end}(1:end-1);       %remove last time step

%sub-segment turns (to better predict beginning and end of saccades)
turnsAll = {};
[turnsAccum,turnsRest] = cutTurns(turns,s);
turnsAll = [turnsAll,turnsAccum];
ntot = 1;
while numel(turnsRest)>0
    [turnsAccum,turnsRest] = cutTurns(turnsRest,s);
    turnsAll = [turnsAll,turnsAccum];
    ntot=ntot+1;
end
turns = turnsAll;



%--------------------- remove problematic saccades -----------------------%
% remove turns where WBF drops below threshold of .001, or where position 
% jumps by more than 4 pixels in a single timestep, or where turn is too
% close to air putt (given by puffInds)
iremove=[];
for i=1:numel(turns)
    if numel(turns{i})==0 || min(wbf(turns{i}))<wbfThresh || any(intersect(puffInds,turns{i}))
        iremove = [iremove,i];
    end
    if max(abs(diff(x(turns{i}))))>4 
        iremove = [iremove,i];
    end
end
turns(iremove)=[];
%-------------------------------------------------------------------------%



%------------------------- find saccades ---------------------------------%
if numel(turns)>0
    for i=1:numel(turns)                                    % iterate for each turn
        ddx(i) = abs(x(turns{i}(end)) - x(turns{i}(1)));    % net change in position (from beginning to end of turn)
        ds(i)  = abs(mean(s(turns{i}))).*ddx(i);            % average WBA times net change in position
        
    end

    M=medcouple(abs(ds)');
    y = quantile(abs(ds),[0.25, 0.5, 0.75]);
    Q3=y(3);
    IQR = iqr(abs(ds));
    thresh = Q3+exp(3*M)*1.5*IQR;   % see https://stats.stackexchange.com/questions/13086/is-there-a-boxplot-variant-for-poisson-distributed-data
    
    ii=find(abs(ds)>thresh);        % these will be defined as saccades
    

    color = zeros(1,120000);
    for i=1:numel(ds)
        if ismember(i,ii)
            color(turns{i}) = abs(ds(i));
        else
            color(turns{i}) = 0;
        end
    end
    
    
else
    ii=[];
    thresh = NaN;

end
%-------------------------------------------------------------------------%


%--------------------- remove problematic saccades -----------------------%
%remove turns in which the direction changes more than once
iremove = [];
for i=1:numel(ii)
    if numel(turns{ii(i)})>0 
        if numel(find(diff(x(turns{ii(i)}))>0))>1 && numel(find(diff(x(turns{ii(i)}))<0))>1
            iremove = [iremove,i];
        end
        
    end
end
ii(iremove)=[];     % these are the indices of the remaining saccades
%-------------------------------------------------------------------------%



%------------------------- extract saccade times -------------------------%

tSaccade = double.empty(0,2);
for count=1:numel(ii)
    iA = turns{ii(count)}(1);
    iB = turns{ii(count)}(end);
    tSaccade = [tSaccade;[iA,iB]];
end
%-------------------------------------------------------------------------%



%---------------------- extract non-saccade times ------------------------%
[~,jj] = sort(tSaccade(:,1));
tSaccade = tSaccade(jj,:);

tnonsaccade = [];
for i=1:size(tSaccade,1)-1
    tnonsaccade = [tnonsaccade;[tSaccade(i,2)+1,tSaccade(i+1,1)-1]];
end
if size(tSaccade,1)>0 && tSaccade(1,1)>1
    tnonsaccade = [[1,tSaccade(1,1)-1];tnonsaccade];
end
if size(tSaccade,1)>0 && tSaccade(end,2)<nT
    tnonsaccade = [tnonsaccade;[tSaccade(end,2)+1,nT]];
end

if size(tSaccade,1)<1
    tnonsaccade = [1,nT];
end

%remove zero-duration fixations
iremove=[];
for i=1:size(tnonsaccade,1)
    inds = tnonsaccade(i,1):tnonsaccade(i,2);
    if numel(inds)==0
        iremove = [iremove,i];
    end
end
tnonsaccade(iremove,:) = [];

%remove regions with airpuffs
ttemp = [];
for i=1:size(tnonsaccade,1)
    inds = tnonsaccade(i,1):tnonsaccade(i,2);
    ipuff = intersect(indsPuff,inds);
    if numel(ipuff)>0
        if (ipuff(1)-1)-tnonsaccade(i,1)>0
            ttemp = [ttemp;[tnonsaccade(i,1),ipuff(1)-1]];
        end
        if tnonsaccade(i,2)-(ipuff(end)+1)>0
            ttemp = [ttemp;[ipuff(end)+1,tnonsaccade(i,2)]];
        end
    else
        ttemp = [ttemp;[tnonsaccade(i,1),tnonsaccade(i,2)]];
    end
end
tnonsaccade = ttemp;


%remove regions with WBF drops (include buffer of 500ms before and after,
%as with airpuffs)
ttemp = [];
for i=1:size(tnonsaccade,1)
    inds = tnonsaccade(i,1):tnonsaccade(i,2);
    
    if numel(find(wbf(inds)<wbfThresh))>30 
        iwbfS = find(wbf(inds)<wbfThresh,1,'first')+inds(1)-1-pmin;
        iwbfF = find(wbf(inds)<wbfThresh,1,'last')+inds(1)-1+pmax;
        if (iwbfS-1)-tnonsaccade(i,1)>0
            ttemp = [ttemp;[tnonsaccade(i,1),iwbfS-1]];
        end
        if tnonsaccade(i,2)-(iwbfF+1)>0
            ttemp = [ttemp;[iwbfF+1,tnonsaccade(i,2)]];
        end
    else
        ttemp = [ttemp;[tnonsaccade(i,1),tnonsaccade(i,2)]];
    end
end
tnonsaccade = ttemp;


%for remaining fixations, compute slope to determine whether fixation or drift:
allSlopes  = [];
allOffsets = [];
allTimes   = [];
tOther     = [];
for i=1:size(tnonsaccade,1)
    inds = tnonsaccade(i,1):tnonsaccade(i,2);
    
    %only keep events longer than 50ms
    if numel(inds)>50

        xin = inds';
        yin = x(inds);
        x0  = xin(1);
        y0  = yin(1);

        P = polyfit(xin-x0,yin-y0,1);
        allSlopes  = [allSlopes; P(1)];
        allOffsets = [allOffsets;P(2)-P(1)*x0+y0];
        allTimes   = [allTimes;[inds(1),inds(end)]]; 
    else
        tOther = [tOther;[inds(1),inds(end)]];
    end
end

threshSlope = .003;

tFixation = [];
sFixation = [];
oFixation = [];

tDrift = [];
sDrift = [];
oDrift = [];
for i=1:numel(allSlopes)
    xtmp = x(allTimes(i,1):allTimes(i,2));
    if var(xtmp)>36 
        tDrift    = [tDrift;   allTimes(i,:)];
        sDrift    = [sDrift;   allSlopes(i) ];
        oDrift    = [oDrift;   allOffsets(i)];
    else
        tFixation = [tFixation;allTimes(i,:)];
        sFixation = [sFixation;allSlopes(i) ];
        oFixation = [oFixation;allOffsets(i)];
    end
end

allIndSac = [];
for i=1:size(tSaccade,1)
    allIndSac = [allIndSac,tSaccade(i,1):tSaccade(i,2)];
end
allIndFix = [];
for i=1:size(tFixation,1)
    allIndFix = [allIndFix,tFixation(i,1):tFixation(i,2)];
end
allIndDrift = [];
for i=1:size(tDrift,1)
    allIndDrift = [allIndDrift,tDrift(i,1):tDrift(i,2)];
end
allIndOther = [];
for i=1:size(tOther,1)
    allIndOther = [allIndOther,tOther(i,1):tOther(i,2)];
end

rsout.fly(fly).trial(trial).saccades.times0    = tSaccade;
rsout.fly(fly).trial(trial).fixations.times0   = tFixation;
rsout.fly(fly).trial(trial).fixations.slopes0  = sFixation;
rsout.fly(fly).trial(trial).fixations.offsets0 = oFixation;
rsout.fly(fly).trial(trial).drift.times0       = tDrift;
rsout.fly(fly).trial(trial).drift.slopes0      = sDrift;
rsout.fly(fly).trial(trial).drift.offsets0     = oDrift;
rsout.fly(fly).trial(trial).other.times0       = tOther;

%----------------------------- plot output -------------------------------%

if plotOutput
    figure;hold on;
    xx = 1:119999;
    ax(1) = subplot(2,2,1);
    scatter(xx,x(1:end-1),20,color(1:end-1),'filled')
    title(['fly: ',num2str(fly),', trial: ',num2str(trial)])
    set(gca,'fontsize',16)
    
    ax(2) = subplot(2,2,3);hold on;set(gca,'fontsize',16)
    kk = find(abs(s)>0);
    cc = color(1:end-1);
    scatter(xx,wbaF(1:end-1),20,color(1:end-1),'filled')
    plot([1,120000],[0,0],'--k');xlim([1,120000])
    
    
    for i=1:numel(ii)
        dr = .1*rand();
        subplot(2,2,3);hold on;plot(turns{ii(i)},dr+3*sign(ds(ii(i)))*ones(1,numel(turns{ii(i)})),'-r','linewidth',2)
    end
    ax(3) = subplot(2,2,2);hold on;set(gca,'fontsize',16)
    t0 = 1:120000;
    x0 = x;
    
    scatter(t0,x0,20,'r','filled')
    scatter(allIndSac,x(allIndSac),20,[0.6902 0.6902 1],'filled')
    for i=1:size(allSlopes,1)
        tfit = allTimes(i,1):allTimes(i,2);
        yfit = allSlopes(i).*tfit+allOffsets(i);
        xtmp = x(allTimes(i,1):allTimes(i,2));
        %if abs(allSlopes(i)) > threshSlope || var(xtmp)>36 %(max(xtmp)-min(xtmp))>28
        if var(xtmp)>36
            scatter(allTimes(i,1):allTimes(i,2),x(allTimes(i,1):allTimes(i,2)),20,[.8,.8,.8],'filled')
            plot(tfit,yfit,'-k','linewidth',2)
        else
            scatter(allTimes(i,1):allTimes(i,2),x(allTimes(i,1):allTimes(i,2)),20,[0 0.7549 0.7941],'filled')
            plot(tfit,yfit,'-','color',[0 0.4549 0.4941],'linewidth',2)
        end
    end

    
    linkaxes(ax,'x')
    
    subplot(2,2,4);hold on;set(gca,'fontsize',16)
    scatter(1:numel(ds),sort(abs(ds)),100,sort(abs(ds)),'filled')
    plot([0,numel(ds)+1],[thresh,thresh],'--k')
    xlim([0,numel(ds)+1])
    title(['number of saccades: ',num2str(numel(ii))])
end
%-------------------------------------------------------------------------%

end


function [turnsAccum,turnsRest] = cutTurns(turns,s)

turnsAccum = {};
turnsRest = {};
m = 1;
p = 1;

%determine whether to break into multiple turns 
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

