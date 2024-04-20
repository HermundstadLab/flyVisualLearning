function plotBehaviorTraces(rsL,rsNL,PIL,PINL,sL,pplot)
% PLOTBEHAVIORTRACES plots the heading trajectories of individual flies and
% shows one trajectory segmented into fixations and saccades
%
% INPUTS:
%   rsNL, rsL: data structures containing the raw heading trajectories
%   sNL, sL:   data structures containining segmented data
%   PINL, PIL:  arrays containing PI scores 
%   pplot: structure containing plotting parameters
%
% See also: PLOTTRAJECTORY, PLOTSACCADES
%

figure('color','w');
    
%laser-trained fly
numFly = 12;
subplot(2,4,1);hold on;plotTrajectory(rsL,numFly,2);set(gca,'fontsize',16);yticks([]);yticklabels({});ylabel('arena position')
title(['laser fly (#',num2str(numFly),'), PI = ', num2str(round(PIL(numFly,2),2))]);xlim([0,120])
subplot(2,4,2);hold on;plotTrajectory(rsL,numFly,3);set(gca,'fontsize',16);yticks([]);yticklabels({});
title(['laser fly (#',num2str(numFly),'), PI = ', num2str(round(PIL(numFly,3),2))]);xlim([0,120])
subplot(2,4,3);hold on;plotTrajectory(rsL,numFly,5);set(gca,'fontsize',16);yticks([]);yticklabels({});
title(['laser fly (#',num2str(numFly),'), PI = ', num2str(round(PIL(numFly,5),2))]);xlim([0,120])
subplot(2,4,4);hold on;plotTrajectory(rsL,numFly,8);set(gca,'fontsize',16);yticks([]);yticklabels({});
title(['laser fly (#',num2str(numFly),'), PI = ', num2str(round(PIL(numFly,5),2))]);xlim([0,120])

%control fly
numFly = 18;
subplot(2,4,5);hold on;plotTrajectory(rsNL,numFly,2);set(gca,'fontsize',16);yticks([]);yticklabels({});ylabel('arena position')
title(['no-laser fly (#',num2str(numFly),'), PI = ', num2str(round(PINL(numFly,2),2))]);xlim([0,120])
xlabel('time (s)')
subplot(2,4,6);hold on;plotTrajectory(rsNL,numFly,3);set(gca,'fontsize',16);yticks([]);yticklabels({});
title(['no-laser fly (#',num2str(numFly),'), PI = ', num2str(round(PINL(numFly,3),2))]);xlim([0,120])
xlabel('time (s)')
subplot(2,4,7);hold on;plotTrajectory(rsNL,numFly,5);set(gca,'fontsize',16);yticks([]);yticklabels({});
title(['no-laser fly (#',num2str(numFly),'), PI = ', num2str(round(PINL(numFly,5),2))]);xlim([0,120])
xlabel('time (s)')
subplot(2,4,8);hold on;plotTrajectory(rsNL,numFly,8);set(gca,'fontsize',16);yticks([]);yticklabels({});
title(['no-laser fly (#',num2str(numFly),'), PI = ', num2str(round(PINL(numFly,8),2))]);xlim([0,120])
xlabel('time (s)')

set(gcf,'Position',[200 200 2000 800])  


numFly   = 12;
numTrial = 2;
figure('color','w');hold on;
plotSaccades(rsL,sL,numFly,numTrial,pplot);
set(gcf,'Position',[200 200 500 350])  
set(gca,'fontsize',16);
xlabel('times (s)')
ylabel('arena position')
