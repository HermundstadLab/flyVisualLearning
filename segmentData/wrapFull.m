function xwrap = wrapFull(x, punishedBar)
% WRAP aligns data to center of danger zone, such that 0 is aligned with
% the center of the danger zone

nx = 96;

if strcmp(punishedBar, 'UpBar')==1
    xwrap = mod(mod(x,nx)-nx/8,nx);
elseif strcmp(punishedBar, 'DnBar')==1
    xwrap = mod(mod(x,nx)-3*nx/8,nx);
elseif strcmp(punishedBar, 'none')==1
    xwrap = mod(x,nx);
else
    error('unrecognized label for punished bar')
end

end