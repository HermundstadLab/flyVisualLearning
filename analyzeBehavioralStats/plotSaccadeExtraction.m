function plotSaccadeExtraction(rsL,sL)
% PLOTSACCADEEXTRACTION illustrates the method for extracting saccades
%
% INPUTS:
%   rsL: data structure containing the raw heading trajectories
%   sL:  data structure containining segmented data
%

%-------------------- plot saccade extraction ------------------------%
extractSaccadesPreprocess(rsL,sL,12,2,'upright','safe',true);
subplot(2,2,1);xlim([56000,67000]);ylim([-30,15]);set(gca,'fontsize',16);xlabel('time (ms)');ylabel('arena position')
subplot(2,2,2);set(gca,'fontsize',16);ylabel('score');xlabel('rank')
subplot(2,2,3);ylim([-2.5,3.5]);set(gca,'fontsize',16);ylabel('\Delta WBA');xlabel('time (ms)')
subplot(2,2,4);xlim([800,910]);ylim([.01,10]);set(gca, 'YScale', 'log','fontsize',16);ylabel('log score');xlabel('rank')

set(gcf,'renderer','Painters','Position',[200 200 1600 1200],'color','w')