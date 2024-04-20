function rC = centerHeadingNeuronFB(theta,heading)
% CENTERHEADINGNEURONFB(THETA, HEADING) generates a sinusoidal tuning curve 
% for a neuron in the FB, defined over THETA and centered at HEADING
%

rC = firingRate(cos(theta - theta(heading)));
end