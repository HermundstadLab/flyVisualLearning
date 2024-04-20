function omega = updateFullPolicy(omega,delta,fk,theta,dx,integ,params,p,policyType)
% UPDATEFULLPOLICY updates the policy of an an artificial agent based on a
% specific set of actions that were taken in the previous timestep
%
% INPUTS:
%   omega: vector of weights
%   delta: reward signal
%   fk: basis functions
%   theta: current heading
%   dx: change in heading
%   integ: an integrator signal used in the drift diffusion process to
%       generate fixations
%   params: a vector of parameters used to specify nonlinear
%       transformations of the weights
%   p: structure containing simulation parameters
%   policyType: string that specifies whether the policy should be trained
%       on an internal representation of heading, or directly on the visual
%       scene. Can take values 'heading' or 'visual'; default is 'heading'
%
% OUTPUTS:
%   omega: updated vector of weights
%

if nargin<9
    policyType = 'heading';
end

nx     = p.nx;
a      = p.aF;
sig    = p.sigF;
delT   = p.delT;
muS    = p.muS;                                      
sigmaS = p.sigmaS;

if strcmp(policyType,'heading')==1
    thetaVec = linspace(0,2*pi,nx+1);
else
    theta = mod(theta,pi);
    thetaVec = linspace(0,pi,nx+1);
end
[~,ind]  = min((theta-thetaVec(1:end-1)).^2);

[nuF,pR,omegaF,omegaR] = convertParams(omega,fk(:,ind),params);

%------------------ FIXATIONS -------------------%

mu  = integ + nuF*delT; 
s   = sig*sqrt(delT);
pF = .5*(1+erf( (a-mu)./sqrt(2.*s.^2) ));

dpF_domega = -delT.*gaussian(a,mu,s).*sigmoidPrime(omegaF*fk(:,ind),params(1),params(2))*fk(:,ind);
dpR_domega = sigmoidPrime(omegaR*fk(:,ind),params(4),params(5))*fk(:,ind);

piFixate   = diracDelta(dx);

%------------------ SACCADES -------------------%

piSaccadeR = lognpdf_branched(dx,muS,sigmaS,'right');
piSaccadeL = lognpdf_branched(dx,muS,sigmaS,'left');

piSaccade  = (1-pR)*piSaccadeL + pR*piSaccadeR;


%------------------ GRADIENTS ------------------%

domegaF = dpF_domega*( piFixate - piSaccade );
domegaR   = (1-pF)*dpR_domega*( piSaccadeR - piSaccadeL );

piEvalGrad = [domegaF;domegaR];    
piEval     = pF*piFixate + (1-pF)*piSaccade;

domega = delta.*piEvalGrad./piEval;
omega  = omega + domega;

end

