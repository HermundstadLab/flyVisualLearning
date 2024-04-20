function j = circIndex(i,nx)
% CIRCINDEX(I, NX) returns the value of the index I over a circular
% interval NX, such that 1<=I<=NX
%
j = mod(i-1,nx)+1;
end