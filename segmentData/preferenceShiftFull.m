function s = preferenceShiftFull(rs, s, p, alignment)
% PREFERENCESHIFT computes the preferred heading from a heading trajectory,
% using the full 360deg of the arena 
% 
% INPUTS:
%   rs: data structure containing the raw heading trajectories
%   s:  data structure containining segmented data
%   p:  parameter vector
%   alignment: a string that can take values 'safe', 'UpBar', 'DnBar', or
%       'none'. Default is 'safe'.
%
% OUTPUT:
%   s: the input data structure, appended with the preference shifts 
%
% See also: PREFERENCESHIFTFULL

if nargin<4
    alignment = 'safe';
end

nx = p.nx;

x0 = nx/4-1;
shifts = 0:(nx-1);
pref = quickWrap(mean(x0)-shifts,nx);
tFlightMin = p.PIindThreshold;

nFlies  = numel(rs.fly);

for i=1:nFlies
    nTrials = numel(rs.fly(i).trial);
    for j=1:nTrials
        
        x    = rs.fly(i).trial(j).x;
        
        if strcmp(alignment,'safe')==1
            xw   = wrapFull(x,rs.fly(i).punishedBar);
        elseif strcmp(alignment,'UpBar')==1
            xw   = wrapFull(x,'DnBar');
        elseif strcmp(alignment,'DnBar')==1
            xw   = wrapFull(x,'UpBar');
        elseif strcmp(alignment,'none')==1
            xw   = wrapFull(x,'none');
        else
            error('unrecognized alignment')
        end
        
        inds = s.fly(i).trial(j).PIinds;
        
        xw = xw(inds);
        
        %use trial iff there is at least 30s of data
        if numel(inds)<tFlightMin
            residencyHist0 = nan(1,nx);
            residencyHist  = residencyHist0;
            PImax = nan;
            dx    = nan;
            pLoc  = nan;

            s.fly(i).trial(j).residencyFull = residencyHist0;
            s.fly(i).trial(j).shiftedFull.PI = PImax;
            s.fly(i).trial(j).shiftedFull.shift = dx;
            s.fly(i).trial(j).shiftedFull.prefLoc = pLoc;
            s.fly(i).trial(j).shiftedFull.residency = residencyHist;

            s.fly(i).trial(j).shiftedFull.saccades.xstartF = [];
            s.fly(i).trial(j).shiftedFull.saccades.xendF   = [];
            s.fly(i).trial(j).shiftedFull.fixations.xavgF  = [];
            
        else
            
            weighting = circshift([zeros(1,nx/4),linspace(0,1,nx/4),fliplr(linspace(0,1,nx/4)),zeros(1,nx/4)],-nx/4);

            PI = zeros(1,numel(shifts));
            f  = zeros(1,numel(shifts));
            for k=1:numel(shifts)
                x = quickWrap(xw + shifts(k),nx);
                PI(k) = (numel([find(x>=nx/8 & x<3*nx/8);find(x>=5*nx/8 & x<7*nx/8)]) ...
                    - numel([find(x<nx/8);find(x>=3*nx/8 & x<5*nx/8);find(x>=7*nx/8)]) )./numel(x);

                c = histcounts(x,0:nx);
                f(k) = sum(weighting.*c);
            end
            %maximize preference within one zone
            [~,imax] = max(f);          
            PImax = PI(imax);
            dx    = shifts(imax);
            pLoc  = pref(imax);

            x = quickWrap(xw + shifts(imax),nx);
            residencyHist = histcounts(x,0:nx);
            residencyHist0 = histcounts(xw,0:nx);

            s.fly(i).trial(j).residencyFull = residencyHist0;
            s.fly(i).trial(j).shiftedFull.PI = PImax;
            s.fly(i).trial(j).shiftedFull.shift = dx;
            s.fly(i).trial(j).shiftedFull.prefLoc = pLoc;
            s.fly(i).trial(j).shiftedFull.residency = residencyHist;

            if isfield(s.fly(i).trial(j).saccades,'dt')
                xstart = quickWrap(s.fly(i).trial(j).saccades.xstartF + dx, nx);
                xend   = quickWrap(s.fly(i).trial(j).saccades.xendF   + dx, nx);

                s.fly(i).trial(j).shiftedFull.saccades.xstartF = xstart;
                s.fly(i).trial(j).shiftedFull.saccades.xendF = xend;
            else
                s.fly(i).trial(j).shiftedFull.saccades.xstartF = [];
                s.fly(i).trial(j).shiftedFull.saccades.xendF = [];
            end


            if isfield(s.fly(i).trial(j).fixations,'dt')
                xavg   = quickWrap(s.fly(i).trial(j).fixations.xavgF + dx, nx);
                s.fly(i).trial(j).shiftedFull.fixations.xavgF = xavg;
            else
                s.fly(i).trial(j).shiftedFull.fixations.xavgF = [];
            end
        end
        
    end
end

end