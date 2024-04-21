%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN SCRIPT                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The following cells are intended to be run in order. Before running:
%
% 1. make sure all folders and subfolders are on the current path
%    analyzeBehavioralStats/:   functions for analyzing base behavioral
%                                   statistics
%    analyzeBumpJumps/:         functions for analyzing calcium imaging
%                                   data 
%    analyzeLearningData/:      functions for analyzing the evolution
%                                   of behavioral data during learning
%    optimizeControlParams/:    functions for optimizing an artificial
%                                   agent with a fully-flexible policy
%    segmentData/:              functions for segmenting behavioral data
%                                   into fixations and saccades
%    simulateCircuitModel/:     functions for simulating an artificial
%                                   agent with a circuit-based model for 
%                                   heading, goal, and policy
%
%    sims/:  folder for storing simulation results
%    fits/:  folder for storing parametric fits to data
%   
% 2. install the MATLAB circular statistics toolbox:
%    https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
%
% 3. set the directory where data is stored:
%
dataDir = 'data/';


%% initiatize parameters for models and plotting

%----------------------- BEHAVIORAL DATA ANALYSIS ------------------------%
p.durThreshold = 0.2;           % minimum duration of fixations in sec
p.velThreshold = 0.1;           % minimum angular velocity of saccades
p.PIindThreshold = 30000;       % minimum amount of time for computing PI scores
p.naiveTrials  = [1,2];         % indices of naive trials
p.trainTrials1 = [3,4];         % indices of first set of training trials
p.probeTrials1 = 5;             % indices of first set of probe trials
p.trainTrials2 = [6,7];         % indices of second set of training trials
p.probeTrials2 = [8,9];         % indices of second set of probe trials
p.nTrials      = 9;             % number of trials
p.nT           = 120000;        % number of timepoints 

%------------------------- FULL MODEL ------------------------------------%
p.tF     = 240000;              % trial duration (in ms)                  
p.nx     = 96;                  % number of pixels in arena
p.nb     = 16;                  % number of basis functions for function approximation
p.aF     = 0.79;                % threshold of DD process for fixations
p.sigF   = 1;                   % spread of DD process for fixations
p.delT   = 0.001;               % timestep (in sec)
p.muS    = 3.89;                % mean param of lognormal distribution over saccades sizes    
p.sigmaS = 0.54;                % spread param of lognormal distribution over saccades sizes 
p.dtS    = 320;                 % duration of saccades

%----------------------- SETPOINT MODEL ----------------------------------%
p.kR     = 1;                   % parameters of sigmoid function for saccades
p.fR     = 0.98;                %                   "           
p.dR     = -0.01;               %                   "       
p.kF     = 1;                   % parameters of sigmoid function for fixations     
p.fF     = 10;                  %                   "
p.dF     = -.01;                %                   "
p.offset = pi/2;                % location of setpoint (in rad)
p.gainF  = .8;                  % gain of fixational setpoint
p.shiftF = .05;                 % shift of fixational setpoint
p.gainS  = .9;                  % gain of saccade setpoint
p.shiftS = 0;                   % shift of saccade setpoint
p.gainB  = .7;                  % gain of bump jump setpoint
p.shiftB = .05;                 % shift of bump jump setpoint
p.shiftB2 = .15;                % alternate shift of bump jump setpoint (for visualization)

%------------------------ CIRCUIT MODEL ----------------------------------%
p.alphaW  = .001;               % learning rate for goal representation
p.alphaH  = .05;                % learning rate for heading representation
p.dtNaive = 30;                 % duration of naive trial (# iterations)
p.dtProbe = 30;                 % duration of probe trial (# iterations)
p.N = 32;                       % discretization of heading space

%-------------------------- PLOTTING -------------------------------------%
pplot.cF = [0,  116,126]./255;
pplot.cS = [125,125,255]./255;
pplot.cB = [174,0  ,66 ]./255;
pplot.cG = [0,  109,243]./255;
pplot.lw = 2;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          EXTRACT DATA                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% load behavioral data from upright arena
% use shorthand names ('rs' for 'results') for notational convenience

rsL    = load([dataDir,'Preproc_AllFly_WT_Laser_FlightBehaviorData.mat'  ]);rsL  = rsL.preProcData_WT_Laser;    %laser-trained flies
rsNL   = load([dataDir,'Preproc_AllFly_WT_NoLaser_FlightBehaviorData.mat']);rsNL = rsNL.preProcData_WT_NoLaser; %no-laser control flies

rsK90  = load([dataDir,'Preproc_AllFly_SS00090_KirSilenced_Laser_FlightBehaviorData.mat']);rsK90 = rsK90.preProcData_KirSS90_Laser; %kir-silenced flies, SS00090
rsK96  = load([dataDir,'Preproc_AllFly_SS00096_KirSilenced_Laser_FlightBehaviorData.mat']);rsK96 = rsK96.preProcData_KirSS96_Laser; %kir-silenced flies, SS00096

rsK90C = load([dataDir,'Preproc_AllFly_SS00090_ParentalControl_Laser_FlightBehaviorData.mat']);rsK90C = rsK90C.preProcData_ParentSS90_Laser;    %parental control flies, SS00090
rsK96C = load([dataDir,'Preproc_AllFly_SS00096_ParentalControl_Laser_FlightBehaviorData.mat']);rsK96C = rsK96C.preProcData_ParentSS96_Laser;    %parental control flies, SS00096

rs2P   = load([dataDir,'Preproc_AllFly_60D05_NoLaser_CaImaging_FlightBehaviorData.mat']);rs2P = rs2P.preProcData_60D05_NoLaser; %two-photon imaging flies

%combine across SS0009X genotypes:
rsKC = rsK90C;
for i=1:numel(rsK96C.fly)
    rsKC.fly(i+numel(rsK90C.fly)) = rsK96C.fly(i);
end

rsK = rsK90;
for i=1:numel(rsK96.fly)
    rsK.fly(i+numel(rsK90.fly)) = rsK96.fly(i);
end


%% load segmented data (or extract saccades/fixations & align based on preference)
% use shorthand names ('s' for 'segmented') for notational convenience

%--------------------- no-laser control flies ----------------------------%
sNL = load([dataDir,'Segmented_AllFly_WT_NoLaser_FlightBehaviorData.mat']);   
sNL = sNL.segData_WT_NoLaser; 

% uncomment to run segmentation from scratch
% sNL = extractSaccadesPreprocess_allFlies(rsNL,'upright');

% tThresh = computeSegmentationThresh(sNL);
% 
% sNL = extractSaccades_allFlies(rsNL,sNL,tThresh);
% sNL = preferenceShift(rsNL,sNL,p);
% sNL = preferenceShiftFull(rsNL,sNL,p);

%----------------------- laser-trained flies -----------------------------%

sL = load([dataDir,'Segmented_AllFly_WT_Laser_FlightBehaviorData.mat']);             
sL = sL.segData_WT_Laser; 

% uncomment to run segmentation from scratch
% sL = extractSaccadesPreprocess_allFlies(rsL,'upright');
% sL = extractSaccades_allFlies(rsL,sL,tThresh);
% sL = preferenceShift(rsL,sL,p);
% sL = preferenceShiftFull(rsL,sL,p);

%-------------------------- SS00090 Kir ----------------------------------%

sK90 = load([dataDir,'Segmented_AllFly_SS00090_KirSilenced_Laser_FlightBehaviorData.mat']);           
sK90 = sK90.segData_KirSS90_Laser; 

% uncomment to run segmentation from scratch
% sK90 = extractSaccadesPreprocess_allFlies(rsK90,'upright');
% sK90 = extractSaccades_allFlies(rsK90,sK90,tThresh);
% sK90 = preferenceShift(rsK90,sK90,p);
% sK90 = preferenceShiftFull(rsK90,sK90,p);

%-------------------------- SS00090 Parent -------------------------------%

sK90C = load([dataDir,'Segmented_AllFly_SS00090_ParentalControl_Laser_FlightBehaviorData.mat']);           
sK90C = sK90C.segData_ParentSS90_Laser; 

% uncomment to run segmentation from scratch
% sK90C = extractSaccadesPreprocess_allFlies(rsK90C,'upright');
% sK90C = extractSaccades_allFlies(rsK90C,sK90C,tThresh);
% sK90C = preferenceShift(rsK90C,sK90C,p);
% sK90C = preferenceShiftFull(rsK90C,sK90C,p);

%---------------------------- SS00096 Kir --------------------------------%

sK96 = load([dataDir,'Segmented_AllFly_SS00096_KirSilenced_Laser_FlightBehaviorData.mat']);           
sK96 = sK96.segData_KirSS96_Laser; 

% uncomment to run segmentation from scratch
% sK96 = extractSaccadesPreprocess_allFlies(rsK96,'upright');
% sK96 = extractSaccades_allFlies(rsK96,sK96,tThresh);
% sK96 = preferenceShift(rsK96,sK96,p);
% sK96 = preferenceShiftFull(rsK96,sK96,p);

%--------------------------- SS00096 Parent ------------------------------%

% uncomment to load segmented data from parental control flies, SS00096 
sK96C = load([dataDir,'Segmented_AllFly_SS00096_ParentalControl_Laser_FlightBehaviorData.mat']);           
sK96C = sK96C.segData_ParentSS96_Laser; 

% uncomment to run segmentation from scratch
% sK96C = extractSaccadesPreprocess_allFlies(rsK96C,'upright');
% sK96C = extractSaccades_allFlies(rsK96C,sK96C,tThresh);
% sK96C = preferenceShift(rsK96C,sK96C,p);
% sK96C = preferenceShiftFull(rsK96C,sK96C,p);


%--------------------------- 60D05 2P imaging --------------------------%

s2P = load([dataDir,'Segmented_AllFly_60D05_NoLaser_CaImaging_FlightBehaviorData.mat']);           
s2P = s2P.segData_60D05_NoLaser; 

% uncomment to run segmentation from scratch
% s2P = extractSaccadesPreprocess_allFlies(rs2P,'2P');
% s2P = extractSaccades_allFlies(rs2P,s2P,tThresh);
% s2P = preferenceShift(rs2P,s2P,p);
% s2P = preferenceShiftFull(rs2P,s2P,p);

%-------------------- combine across SS0009X genotypes -------------------%
sK = sK90;
for i=1:numel(sK96.fly)
    sK.fly(i+numel(sK90.fly)) = sK96.fly(i);
end

sKC = sK90C;
for i=1:numel(sK96C.fly)
    sKC.fly(i+numel(sK90C.fly)) = sK96C.fly(i);
end

% uncomment to save segmented data:
% segData_WT_Laser         = sL;
% segData_WT_NoLaser       = sNL;
% segData_KirSS90_Laser    = sK90;
% segData_KirSS96_Laser    = sK96;
% segData_ParentSS90_Laser = sK90C;
% segData_ParentSS96_Laser = sK96C;
% segData_60D05_NoLaser  = s2P;
% save([dataDir,'Segmented_AllFly_WT_Laser_FlightBehaviorData.mat'],'segData_WT_Laser');
% save([dataDir,'Segmented_AllFly_WT_NoLaser_FlightBehaviorData.mat'],'segData_WT_NoLaser');
% save([dataDir,'Segmented_AllFly_SS00090_KirSilenced_Laser_FlightBehaviorData.mat'],'segData_KirSS90_Laser');
% save([dataDir,'Segmented_AllFly_SS00090_ParentalControl_Laser_FlightBehaviorData.mat'],'segData_ParentSS90_Laser');
% save([dataDir,'Segmented_AllFly_SS00096_KirSilenced_Laser_FlightBehaviorData.mat'],'segData_KirSS96_Laser');
% save([dataDir,'Segmented_AllFly_SS00096_ParentalControl_Laser_FlightBehaviorData.mat'],'segData_ParentSS96_Laser');
% save([dataDir,'Segmented_AllFly_60D05_NoLaser_CaImaging_FlightBehaviorData.mat'],'segData_60D05_NoLaser');

%% process 2P imaging data

ttmp = load([dataDir,'ProcJumpData_AllFly_60D05_NoLaser_CaImaging_FlightBehaviorData.mat']); 
tSeriesProcAll = ttmp.tSeriesProcAll; 

% uncomment to process and save bump jump data:
% tSeriesProcAll = ProcessEBCaImagingAndFlightBehaviorData(rs2P, s2P);
% save([dataDir,'ProcJumpData_AllFly_60D05_NoLaser_CaImaging_FlightBehaviorData.mat'], 'tSeriesProcAll','-v7.3');


%% fit parametric functions to distributions of fixations and saccades
% fitting results are stored in subfolder /fits/
fitFixDists(sL,sNL);            %check that IG is best-fitting distribution
fitFixationSaccades(sL,sNL,p);  %fit IG distribution to data


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GENERATE PLOTS                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%---------------------------- FIG 1 ------------------------------------%
%
%% plot increase in behavioral preference over time (Fig 1F & SI Fig S2A,D-F)
plotPIscores(rsL,rsNL,sL,sNL,'laser','trialPref',p);        %compute PI score wrt preference in each trial (Fig 1F)


%% plot standard PI scores (Fig 1G-H)
[PIL,PINL] = plotPIscores(rsL,rsNL,sL,sNL,'laser','both',p);    %Fig 1G
[PIKC,PIK] = plotPIscores(rsKC,rsK,sKC,sK,'kir','both',p);      %Fig 1H


%% plot PI scores split out by pattern (Fig 1J-K )
plotPIscores(rsKC,rsK,sKC,sK,'kir','UpBar',p);  %Fig 1J
plotPIscores(rsKC,rsK,sKC,sK,'kir','DnBar',p);  %Fig 1K


%% plot behavioral trajectories (Fig 1D-E, Fig 3A)
plotBehaviorTraces(rsL,rsNL,PIL,PINL,sL,pplot);


%%
%------------- SI FIGS S1, S2, S3, & S4 (related to FIG 1) ---------------%
%
%% SI Fig S1: plot heading preferences
plotHeadingPreferences(sL,sNL,rsL,rsNL,p);


%% SI Fig S2: plot PI scores

%plot distribution of airpuffs 
plotAirpuffs(rsL,sL,p);                                     %SI Fig S2C   

%plot increase in behavioral preference over time
plotPIscores(rsL,rsNL,sL,sNL,'laser','finalPref',p);        %SI Fig S2A (compute PI score wrt final preference in trial 9)
plotPIscores(rsKC,rsK,sKC,sK,'kir','trialPref',p);          %SI Fig S2D 
plotPIscores(rsK90C,rsK90,sK90C,sK90,'kir','trialPref',p);  %SI Fig S2E  
plotPIscores(rsK96C,rsK96,sK96C,sK96,'kir','trialPref',p);  %SI Fig S2F

%plot standard PI scores
plotPIscores(rsL,rsNL,sL,sNL,'laser','both',p,1);           %SI Fig S2B (exclude flies that received airpuffs) 
plotPIscores(rsK90C,rsK90,sK90C,sK90,'kir','both',p);       %SI Fig S2G 
plotPIscores(rsK96C,rsK96,sK96C,sK96,'kir','both',p);       %SI Fig S2J 

%plot PI scores split out by pattern
plotPIscores(rsK90C,rsK90,sK90C,sK90,'kir','DnBar',p);      %SI Fig S2H
plotPIscores(rsK90C,rsK90,sK90C,sK90,'kir','UpBar',p);      %SI Fig S2I 

plotPIscores(rsK96C,rsK96,sK96C,sK96,'kir','DnBar',p);      %SI Fig S2K
plotPIscores(rsK96C,rsK96,sK96C,sK96,'kir','UpBar',p);      %SI Fig S2L 

%% SI Fig S3: perform ANOVA
AOV = performANOVA(rsL,rsNL,rsK90C,rsK90,rsK96C,rsK96,sL,sNL,sK90C,sK90,sK96C,sK96,p);


%% SI Fig S4: 
plotStripeFixation(rsK90C,rsK90,rsK96C,rsK96,rsL,rsNL);


%%
%---------------------------- FIG 2 --------------------------------------%
%
%% plot final compass mapping (Fig 2B-D)
illustrateVisualMapping(p,pplot,1,{'asym',24,.5})         %asymmetric scene
illustrateVisualMapping(p,pplot,1,{'symNoJump', 24,.5})   %symmetric scene, assuming no jumps


%% plot final heading weights (Fig 2G)
illustrateVisualMapping(p,pplot,2)    %simulate heading trajectory with jumps


%% plot heading trajectory for fly 20 (FIG 2H)
plotHeadingTrajectory(tSeriesProcAll,20)


%% plot offsets for all flies (FIG 2I, left panel)
plotOffsetsAcrossFlies(tSeriesProcAll,1:9);


%% plot summary of offsets (FIG 2I, right panels)
[~,~,nOffsets,~] = plotOffsetSummary(tSeriesProcAll,1:9);


%% plot firing rate map (FIG 2J)
plotDoublyPeakedTuningCurves(tSeriesProcAll,20,1:9)


%%
%---------------------------- FIG 3 --------------------------------------%
%
%% plot segmented behavioral trajectory (Fig 3A, matched to Fig 1D-E)
plotBehaviorTraces(rsL,rsNL,PIL,PINL,sL,pplot);


%% plot optimal duration of fixations and direction of saccades as a function of heading (Fig 3D,E)
plotOptimalPolicy(p,pplot,'heading');


%% plot heading-dependent averages, aligned to preference (Fig 3F)
plotBehaviorAvgs(sL, 'full_flat','pref',p,pplot,5,p.naiveTrials);     %laser-trained, naive trials
plotBehaviorAvgs(sNL,'full_flat','pref',p,pplot,5,p.naiveTrials);     %no-laser controls, naive trials

probeInds = [p.trainTrials1,p.probeTrials1,p.trainTrials2,p.probeTrials2];
plotBehaviorAvgs(sL, 'full_flat','pref',p,pplot,5,probeInds);         %laser-trained, train/probe trials
plotBehaviorAvgs(sNL,'full_flat','pref',p,pplot,5,probeInds);         %no-laser controls, train/probe trials


%% plot evolution of the goal heading (Fig 3G)
pvals = plotLearningEvolution(1,sL,rsL,sNL,rsNL,[],[],p);


%%
%------------------- SI FIGS S5-S6 (related to FIG 3) ------------------%
%
%% SI Fig S5: segment behavior

% extract saccades
plotSaccadeExtraction(rsL,sL)   % SI Fig S5A

% plot joint distribution of fixations and saccades 
plotJointVelDurDist(sL,sNL,p);  % SI Fig S5B

% plot empirical distributions of fixations and saccades
plotEmpiricalPDFs(sL,sNL,p);    % SI Fig S5C-F

% generate saccades-size parameters used in model
[medDur,locP,scaleP,locD,scaleD] = approxSaccadeParams(sL,sNL);

% plot inverse Gaussian fits to duration distribtions
plotFixSacFits(sL,p);           % SI Fig S5H-I


%% SI Fig S6: plot heading-dependent averages

probeInds = [p.trainTrials1,p.probeTrials1,p.trainTrials2,p.probeTrials2];

% plot heading dependent averaged, centered on safety
plotBehaviorAvgs(sL, 'half','safe',p,pplot,5,p.naiveTrials,'saccadeSize');     % SI Fig S6A: laser-trained, naive trials
plotBehaviorAvgs(sNL,'half','safe',p,pplot,5,p.naiveTrials,'saccadeSize');     % SI Fig S6A: no-laser controls, naive trials

plotBehaviorAvgs(sL, 'half','safe',p,pplot,5,probeInds,'saccadeSize');         % SI Fig S6A: laser-trained, train/probe trials
plotBehaviorAvgs(sNL,'half','safe',p,pplot,5,probeInds,'saccadeSize');         % SI Fig S6A: no-laser controls, train/probe trials

% plot collapsed heading-dependent averages, centered on preference
plotBehaviorAvgs(sL, 'half','pref',p,pplot,5,p.naiveTrials,'saccadeSize');     % SI Fig S6B: laser-trained, naive trials
plotBehaviorAvgs(sNL,'half','pref',p,pplot,5,p.naiveTrials,'saccadeSize');     % SI Fig S6B: no-laser controls, naive trials

plotBehaviorAvgs(sL, 'half','pref',p,pplot,5,probeInds,'saccadeSize');         % SI Fig S6B: laser-trained, train/probe trials
plotBehaviorAvgs(sNL,'half','pref',p,pplot,5,probeInds,'saccadeSize');         % SI Fig S6B: no-laser controls, train/probe trials

%to generate averages of fixation WBF, saccade velocity or saccade WBF,
%repeat the above with one of the following final inputs (instead of
%'saccadeSize'): 'fixationFrequency, 'saccadeVelocity', 'saccadeFrequency'


%%
%---------------------------- FIG 4 --------------------------------------%
%
%% plot illustration of multiplicative operation (Fig 4C-D)
illustratePolicy(p,pplot,1);


%% plot illustration of phase shift + multiplicative operation (Fig 4E)
illustratePolicy(p,pplot,2)


%% plot illustration of PFL activity (Fig 4F)
illustratePolicy(p,pplot,3)


%% illustration of goal update  (Fig 4G)
illustratePolicy(p,pplot,4);


%% plot evolution of goal heading for fixed compass heading  (Fig 4H)
illustratePolicy(p,pplot,5);


%%
%---------------------------- FIG 5 --------------------------------------%
%
%% load data
load('sims/PIsweep_sym_probe_uniform_asym_init.mat')
%sweepPIscores(p,{'sym','probe','uniform','asym_init'});    % uncomment to generate simulations with fixed goal (uses parallel pool)


%% plot heading weights over time (Fig 5C), conditional probability of a jump (Fig 5D, left panel), and final goal/map strengths (Fig 5E)
plotSims_VisualMapping(p,headGains,headingJumps,headingVisits,PIscores,gain,1,{'sym','probe','asym_init'});
plotSims_VisualMapping(p,headGains,headingJumps,headingVisits,PIscores,gain,1,{'sym','probe','rand_init'});


%% plot best fitting cosine curve for bump jump probability (FIG 5D, right panel)
plotProbBumpJump(tSeriesProcAll,nOffsets,p);


%%
%-------------------- SI FIG S8 (related to FIG 5) ---------------------%
%
%% SI Fig S8: different mapping conditions 

%load data: asymmetric initial map
load('sims/PIsweep_sym_probe_uniform_asym_init.mat');
%sweepPIscores(p,{'sym','probe','uniform','asym_init'});    % uncomment to generate simulations with fixed goal (uses parallel pool)

plotSims_VisualMapping(p,headGains,headingJumps,headingVisits,PIscores,gain,1,{'sym','probe','asym_init'});   % SI Fig S8A-C

%load data: random initial map
load('sims/PIsweep_sym_probe_uniform_rand_init.mat');
%sweepPIscores(p,{'sym','probe','uniform','rand_init'});    % uncomment to generate simulations with fixed goal (uses parallel pool)

plotSims_VisualMapping(p,headGains,headingJumps,headingVisits,PIscores,gain,1,{'sym','probe','randinit'});    % SI Fig S8D-F
    

%%
%---------------------------- FIG 6 --------------------------------------%
%
%% load data for uniform initial goal location
load('sims/PIsweep_sym_train_uniform_asym_init.mat')
%sweepPIscores(p,{'sym','train','uniform','asym_init'});      % generate simulations with evolving goal (uses parallel pool)


%% generate vector maps (Fig 6D)
plotSims_GoalLearning_VisualMapping(1,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);


%% plot individual fly trajectories (Fig 6E-F)
plotSims_GoalLearning_VisualMapping(2,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);


%% plot individual variability across learning simulations (Fig 6G-H)
plotSims_GoalLearning_VisualMapping(3,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);


%% plot variance in goal updates (Fig 6G and SI Fig S9)
plotSims_GoalLearning_VisualMapping(4,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);


%% plot fraction of flies with different initial goal headings (Fig 6G)
plotSims_GoalLearning_VisualMapping(5,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);


%% load data for random initial goal location
load('sims/PIsweep_sym_train_random_asym_init.mat')
%sweepPIscores(p,{'sym','train','random','asym_init'});      % generate simulations with evolving goal (uses parallel pool)


%% plot evolution over time, model (Fig 6I, left column)
plotSims_GoalLearning_VisualMapping(6,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);


%% plot evolution over time, data (Fig 6I, right column)
pvals = plotLearningEvolution(2,sL,rsL,sNL,rsNL,[],[],p);


%%
%-------------------- SI FIGS S9, S10 (related to FIG 6) ---------------%
%
%% SI Fig S9: plot variance in goal updates

%load data
load('sims/PIsweep_sym_train_uniform_asym_init.mat')
%sweepPIscores(p,{'sym','train','random','asym_init'});      % uncomment to generate simulations with evolving goal (uses parallel pool)

plotSims_GoalLearning_VisualMapping(4,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);


%% SI Fig S10: plot evolution over time, grouped in different ways

%load data: random initial conditions
load('sims/PIsweep_sym_train_random_rand_init.mat')
%sweepPIscores(p,{'sym','train','random','rand_init'});      % uncomment to generate simulations with evolving goal (uses parallel pool)

%plot CDF of initial conditions
plotSims_GoalLearning_VisualMapping(7,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);  % SI Fig S10A: model
plotLearningEvolution(4,sL,rsL,sNL,rsNL,sKC,rsKC,p,0);                             % SI Fig S10A: data

%plot evolution over time, model
plotSims_GoalLearning_VisualMapping(6,headLocs,headGains,goalLocs,goalGains,distDanger,PIscores,angdiff,gain);  % SI Fig S10B: model

%plot evolution over time, data
pvals_twoGroups   = plotLearningEvolution(2,sL,rsL,sNL,rsNL,sKC,rsKC,p,1);         % SI Fig S10C: data
pvals_threeGroups = plotLearningEvolution(3,sL,rsL,sNL,rsNL,sKC,rsKC,p,1);         % SI Fig S10D: data     

