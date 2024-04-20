function [dx,dt,integ] = sampleFullPolicy(omega,theta,fk,integ,params,p)
% SAMPLEFULLPOLICY samples a change in orientation and over a change in 
% time to use in training an artificial agent.
%
% INPUTS:
%   omega: vector of weights
%   theta: current heading
%   fk: basis functions
%   integ: an integrator signal used in the drift diffusion process to
%       generate fixations
%   params: a vector of parameters used to specify nonlinear
%       transformations of the weights
%   p: structure containing simulation parameters
%
% OUTPUTS: 
%   dx: change in orientation
%   dt: change in time
%   integ: updated integratory variable
%

nx       = p.nx;
thetaVec = linspace(0,2*pi,nx+1);
[~,ind]  = min((theta-thetaVec(1:end-1)).^2);

%update integrator
[nu,pR] = convertParams(omega,fk(:,ind),params);
integ   = integ + normrnd(nu*p.delT,p.sigF*sqrt(p.delT));

if integ>p.aF
    %--------- SACCADE ----------%
    if flipCoin(pR)
        dir = 1;
    else
        dir = -1;
    end
    
    dx = dir.*lognrnd(p.muS,p.sigmaS)*96./360;
    dt = p.dtS;
    integ = 0;

else 
    %--------- FIXATE ----------%
    dx = 0;
    dt = 1;
end

dx = round(dx);
dt = round(dt);

end