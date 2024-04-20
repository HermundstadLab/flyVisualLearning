function x = scale(x,dx)
% SCALE(X,DX) rescales an input profile X (spanning 0 to 1) to span DX to
% (1-DX)
%
x = (1-2*dx)*x+dx;
end