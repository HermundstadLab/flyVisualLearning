function x = normalize(x)
% NORMALIZE(X) normalizes the activity profile X to be between 0 and 1
%
x = (x-min(x))./(max(x)-min(x));
end