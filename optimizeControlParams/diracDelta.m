function p = diracDelta(x0)
% DIRACDELTA(X0)=delta(X0-0) generates the output of dirac delta function 
% centered at zero 
%

if abs(x0)<1e-10
    p = 1;
else
    p = 0;
end