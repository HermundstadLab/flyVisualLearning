function fk = getBasisFunctions(nb,k,nx)
% GETBASISFUNCTIONS generates a set of von Mises basis functions that evenly 
% tile a circular interval
%
% INPUTS:
%   nb: number of basis functions
%   k: concentration of von Mises profile
%   nx: number of bins to use to discretize of input space
%


x=linspace(-pi/2,pi/2,nx);
centers = linspace(-pi/2,pi/2,nb+1);
f0 = vonMises(x,0,k);
[~,i0] = min(abs(x-0));

fk=zeros(nb,numel(x)-1);
for i=1:nb
    [~,j] = min(abs(x-centers(i)));
    if i0>=j
        fk(i,:) = [f0(i0-j+1:end-1),f0(1:i0-j)];
    else
        fk(i,:) = [f0(end-(j-i0):end-1),f0(1:end-(j-i0)-1)];
    end
end
end