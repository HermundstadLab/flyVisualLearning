function s = preferenceShift(rs, s, p, alignment)
% PREFERENCESHIFT computes the preferred heading from a heading trajectory,
% collapsed onto 180deg of the arena using the two-fold symmetry of the
% visual environment
% 
% INPUTS:
%   rs: data structure containing the raw heading trajectories
%   s:  data structure containining segmented data
%   alignment: a string that can take values 'safe', 'UpBar', 'DnBar', or
%       'none'. Default is 'safe'.
%
% OUTPUT:
%   s: the input data structure, appended with the preference shifts 
%
% See also: PREFERENCESHIFTFULL
%

if nargin<4
    alignment = 'safe';
end

nx = p.nx/2;

shifts = 0:(nx-1);
x0 = (nx/2-1);
pref = quickWrap(mean(x0)-shifts, nx);

tFlightMin = p.PIindThreshold;

nFlies  = numel(rs.fly);

for i=1:nFlies
    nTrials = numel(rs.fly(i).trial);
    for j=1:nTrials
        x    = rs.fly(i).trial(j).x;
        
        if strcmp(alignment,'safe')==1
            xw   = wrap(x,rs.fly(i).punishedBar);
        elseif strcmp(alignment,'UpBar')==1
            xw   = wrap(x,'DnBar');
        elseif strcmp(alignment,'DnBar')==1
            xw   = wrap(x,'UpBar');
        elseif strcmp(alignment,'none')==1
            xw   = wrap(x,'none');
        else
            error('unrecognized alignment')
        end
        
        inds = s.fly(i).trial(j).PIinds; 
        
        xw = xw(inds);
        
        if numel(inds)<tFlightMin
            residencyHist  = nan(1,nx);
            residencyHist0 = residencyHist;
            PImax = nan;
            dx    = nan;
            pLoc  = nan;
            
            s.fly(i).trial(j).residency = residencyHist0;
            s.fly(i).trial(j).shifted.PI = PImax;
            s.fly(i).trial(j).shifted.shift = dx;
            s.fly(i).trial(j).shifted.prefLoc = pLoc;
            s.fly(i).trial(j).shifted.residency = residencyHist;

            s.fly(i).trial(j).shifted.saccades.xstart = [];
            s.fly(i).trial(j).shifted.saccades.xend   = [];
            s.fly(i).trial(j).shifted.fixations.xavg  = [];
            
        else
        

            weighting = [linspace(0,1,nx/2),fliplr(linspace(0,1,nx/2))];
            
            PI = zeros(1,numel(shifts));
            f  = zeros(1,numel(shifts));
            for k=1:numel(shifts)
                x = quickWrap(xw + shifts(k), nx);
                PI(k) = (numel([find(x>=nx/4 & x<3*nx/4);find(x>=5*nx/4 & x<7*nx/4)]) ...
                    - numel([find(x<nx/4);find(x>=3*nx/4 & x<5*nx/4);find(x>=7*nx/4)]) )./numel(x);
                
                c = histcounts(x,0:nx);
                f(k) = sum(weighting.*c);
            end
            %maximize preference within one zone
            [~,imax] = max(f);         
            PImax = PI(imax);
            dx    = shifts(imax);
            pLoc  = pref(imax);

            x = quickWrap(xw + shifts(imax), nx);
            residencyHist = histcounts(x,0:nx);


            residencyHist0 = histcounts(xw,0:nx);
            
            s.fly(i).trial(j).residency = residencyHist0;
            s.fly(i).trial(j).shifted.PI = PImax;
            s.fly(i).trial(j).shifted.shift = dx;
            s.fly(i).trial(j).shifted.prefLoc = pLoc;
            s.fly(i).trial(j).shifted.residency = residencyHist;

            if isfield(s.fly(i).trial(j).saccades,'dt')
                xstart = quickWrap(s.fly(i).trial(j).saccades.xstart + dx, nx);
                xend   = quickWrap(s.fly(i).trial(j).saccades.xend   + dx, nx);

                s.fly(i).trial(j).shifted.saccades.xstart = xstart;
                s.fly(i).trial(j).shifted.saccades.xend = xend;
            else
                s.fly(i).trial(j).shifted.saccades.xstart = [];
                s.fly(i).trial(j).shifted.saccades.xend = [];
            end


            if isfield(s.fly(i).trial(j).fixations,'dt')
                xavg   = quickWrap(s.fly(i).trial(j).fixations.xavg + dx, nx);
                s.fly(i).trial(j).shifted.fixations.xavg = xavg;
            else
                s.fly(i).trial(j).shifted.fixations.xavg = [];
            end
        end
        
        
    end
end

end



