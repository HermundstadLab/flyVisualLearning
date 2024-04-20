function [omega,fk,params,xAll,nuAll,prAll] = trainFullPolicy(p,policyType)
% TRAINFULLPOLICY trains an artificial agent using a learnable set of
% weights onto a set of basis functions. This set of weights is then used
% to specify the drift rate of fixations and the directional bias of
% saccades at different headings.
%
% INPUTS:
%   p: structure containing simulation parameters
%   policyType: string that specifies whether the policy should be trained
%       on an internal representation of heading, or directly on the visual
%       scene. Can take values 'heading' or 'visual'; default is 'heading'
%
% OUTPUTS:
%   omega: vector of weights
%   fk: basis functions
%   params: a vector of parameters used to specify nonlinear
%       transformations of the weights
%   xAll: the trajectory of the agent over time
%   nuAll: the drift rates of the agent over time
%   prAll: the probability of rightward saccades over time
%

if nargin<2
    policyType = 'heading';
end

%parameters for basis functions
nb = p.nb;      % number of basis functions (8 or 16 in fly EB)
k  = nb/2;      % width of basis functions
nx = p.nx;      % discretization of stimulus space (corresponds to number of pixels in visual arena)
fk = getBasisFunctions(nb,k,nx+1);

%simulation parameters
tF    = p.tF/4;                                 % total time  
t     = 1;                                      % initial time
x     = randperm(nx,1);                         % linear position in arena
theta = mod(2*pi*(x-1)./nx,2*pi);

%parameters for saccades
kR   = p.kR;       % slope of nonlinearity
fR   = p.fR;       % amplitude of nonlinearity (corresponds to Pmax - Pmin = 1)
dR   = p.dR;       % offset of nonlinearity (corresponds to Pmax = 1)

%policy for fixations: choose params s.t. nuMin = .1, nuMax = .5,                
kF   = p.kF;       % slope of nonlinearity
fF   = p.fF;       % amplitude of nonlinearity (corresponds to nuMax-nuMin = .4)
dF   = p.dF;       % offset of nonlinearity (corresponds to nuMin = .1)

%policy parameters
omegaF = .1;
omegaS = 0;
omega  = [omegaF*ones(nb,1);...           % learnable fixation params
    omegaS*ones(nb,1)];                   % learnable saccade params 

params = [kF,fF,dF,kR,fR,dR];   
integ  = 0;

[nuF,pR,~,~] = convertParams(omega,fk,params);
xAll = x;
nuAll = nuF;
prAll = pR;


while t<tF

    %policy is tethered to bump position, thetaBump
    [dx,dt,integnew] = sampleFullPolicy(omega,theta,fk,integ,params,p);

    %action is taken based on virtual position, xVirtual
    tnew     = t+dt; 
    xnew     = x+dx;
    thetanew = mod(2*pi*(xnew-1)./nx,2*pi);

    %compute RPE
    delta = getDelta(thetanew,policyType);

    %update parameters
    omega = updateFullPolicy(omega,delta,fk,theta,dx,integ,params,p,policyType);

    [nuF,pR,~,~] = convertParams(omega,fk,params);
    if dx==0
        xAll  = [xAll,x+normrnd(0,.5)];
        nuAll = [nuAll,nuF];
        prAll = [prAll,pR];
    else
        xAll  = [xAll,linspace(x,x+dx,round(dt))];
        nuAll = [nuAll,repmat(nuF,[1,round(dt)])];
        prAll = [prAll,repmat(pR,[1,round(dt)])];
    end
    
    %update states
    t     = tnew; 
    x     = xnew;      
    theta = thetanew;
    integ = integnew;

end
nuAll = reshape(nuAll',[nx,numel(nuAll)./nx])';
prAll = reshape(prAll',[nx,numel(prAll)./nx])';

    
end