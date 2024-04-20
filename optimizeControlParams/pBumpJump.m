function pJump = pBumpJump(theta,offset,gain,shift)
% PBUMPJUMP(THETA, OFFSET, GAIN, SHIFT)= (GAIN/2)*(cos(THETA-OFFSET-pi)+1)+SHIFT 
% specifies the probability of a bump jump at an orientation THETA. The 
% OFFSET, GAIN, and SHIFT are used to determine location and scale of the 
% sinusoidal function that specifies this probability.

pJump = gain./2.*(cos(theta-offset-pi)+1) + shift;
