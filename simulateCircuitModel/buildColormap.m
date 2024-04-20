function cmap = buildColormap(c1,c2,np)
%BUILDCOLORMAP(C1,C2,NP) uses NP points to interpolate between colors C1 
% and C2 to build a custom colormap. This assumes C2 is lighter (larger)
% than C1
%

cmap = [];
for i=1:np
    cmap = [cmap;c1+i*(c2-c1)./np];
end
cmap = flipud(cmap);

end