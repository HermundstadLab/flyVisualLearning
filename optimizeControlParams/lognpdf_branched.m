function p = lognpdf_branched(x,mu,sig,branch)
% LOGNPDF_BRANCHED generates the probability of an input, as specified by a 
% lognormal distribution defined on the negative or positive axis. 
%
% INPUTS:
%   x: input
%   mu, sig: parameters of lognormal distribution
%   branch: string that specifies whether inputs are distributed over the
%       positive or negative axis
%


if strcmp(branch,'left')
    p = sign(x) * (sign(x)-1)/2 * lognpdf(-x,mu,sig);

elseif strcmp(branch,'right')
    p = (sign(x)+1)/2 * lognpdf(x,mu,sig); 
end

end