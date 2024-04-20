function tSeriesProcAll = ProcessEBCaImagingAndFlightBehaviorData(results2PBehavior, s2P)
% Process two-photon Ca imaging data from EPG neurons in the 
% Ellipsoid Body (EB) of tethered flying flies. Combine with behavioral
% data.
%
% Vivek Jayaraman

% Reformat the EB (EPG) imaging data and setpoint behavioral data
procData = ReformatEPGAndBehaviorData(results2PBehavior, s2P);
tSeriesProcAll = CompileEPGBumpJumpDataAllFlies(procData);

end % of function
