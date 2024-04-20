function xw = quickWrap(x,nx)
% QUICKWRAP(X,NX) wraps the data in X based on the interval NX, such that
% 0<=X<NX

xw = x;
xw(x>(nx-1)) = x(x>(nx-1))-nx;
xw(x<0) = x(x<0)+nx;
end