function rC = centerHeadingNeuronEB(theta,heading)
% CENTERHEADINGNEURONEB(THETA, HEADING) generates a von Mises tuning curve 
% for a neuron in the EB, defined over THETA and centered at HEADING
%

rC = vonMises(theta,theta(heading),pi);
rC = scale(normalize(rC),.05);
end