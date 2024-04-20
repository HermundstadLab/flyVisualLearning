function coin = flipCoin(p)
% FLIPCOIN(P) samples from a Bernoulli process with probability P
%

u = rand();
if u<p
    coin = true;
else
    coin = false;
end
end