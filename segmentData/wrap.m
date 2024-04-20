function xwrap = wrap(x,punishedBar)
% WRAP aligns collapsed data to center of danger zone, such that 0 is 
% aligned with the center of the danger zone

nx = 48;

if strcmp(punishedBar, 'UpBar')==1
    xwrap = mod(mod(x,nx)-nx/4,nx);
elseif strcmp(punishedBar, 'DnBar')==1
    xwrap = mod(mod(x,nx)-3*nx/4,nx);
elseif strcmp(punishedBar, 'none')==1
    xwrap = mod(x,nx);
else
    error('unrecognized label for punished bar')
end

end