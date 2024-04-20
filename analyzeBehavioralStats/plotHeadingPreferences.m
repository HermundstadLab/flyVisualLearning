function plotHeadingPreferences(sL,sNL,rsL,rsNL,p)
% PLOTHEADINGPREFERENCES plots the heading preferences of individual flies
%
% INPUTS:
%   rsL, rsNL, rsKC: data structures containing the raw heading trajectories
%   sL, sNL:  data structures containining segmented data
%   p: structure containing parameter values
%


nx = p.nx/2;
for i=1:numel(sL.fly)
    for j=1:numel(sL.fly(i).trial)
        x  = rsL.fly(i).trial(j).x;
        xw = wrap(x,rsL.fly(i).punishedBar);

        inds = sL.fly(i).trial(j).PIinds; 
        xw = xw(inds);
        x  = quickWrap(xw + sL.fly(i).trial(j).shifted.shift,nx);
            
        rL(i,j) = (numel(find(x>=nx/4 & x<3*nx/4))-numel([find(x<nx/4);find(x>=3*nx/4)]))./numel(x);    
        prefL(i,j) = sL.fly(i).trial(j).shiftedFull.prefLoc*pi./48;
        nindsL(i,j) = numel(sL.fly(i).trial(j).PIinds);
    end
end

for i=1:numel(sNL.fly)
    for j=1:numel(sNL.fly(i).trial)
        x  = rsNL.fly(i).trial(j).x;
        xw = wrap(x,rsNL.fly(i).punishedBar);

        inds = sNL.fly(i).trial(j).PIinds; 
        xw = xw(inds);
        x  = quickWrap(xw + sNL.fly(i).trial(j).shifted.shift,nx);
            
        rNL(i,j) = (numel(find(x>=nx/4 & x<3*nx/4))-numel([find(x<nx/4);find(x>=3*nx/4)]))./numel(x);
        prefNL(i,j) = sNL.fly(i).trial(j).shiftedFull.prefLoc*pi./48;
        nindsNL(i,j) = numel(sNL.fly(i).trial(j).PIinds);
    end
end

cmap = parula(9);

trial = 2;
figure;set(gcf,'Position',[200 200 1400 600],'color','w')
subplot(2,4,1);polarplot([0,prefL(1,trial)],[0,rL(1,trial)],'color',cmap(trial-1,:));
hold on;
for i=2:numel(rsL.fly)
    polarplot([0,prefL(i,trial)],[0,rL(i,trial)],'color',cmap(trial-1,:))
    rlim([0,1])
    rticks([.5,1])
    set(gca,'fontsize',16)
end
title('laser-trained, trial 2')

subplot(2,4,2);polarplot([0,prefNL(1,trial)],[0,rNL(1,trial)],'color',cmap(trial-1,:));
hold on;
for i=2:numel(rsNL.fly)
    polarplot([0,prefNL(i,trial)],[0,rNL(i,trial)],'color',cmap(trial-1,:))
    rlim([0,1])
    rticks([.5,1])
    set(gca,'fontsize',16)
end
title('no-laser controls, trial 2')

trial = 7;
subplot(2,4,3);polarplot([0,prefNL(1,trial)],[0,rNL(1,trial)],'color',cmap(trial-1,:));
hold on;
for i=2:40
    polarplot([0,prefNL(i,trial)],[0,rNL(i,trial)],'color',cmap(trial-1,:))
    rlim([0,1])
    rticks([.5,1])
    set(gca,'fontsize',16)
end
title('no-laser controls, trial 7')

for trial = 2:9
    subplot(2,4,4);polarplot([0,prefNL(1,trial)],[0,rNL(1,trial)],'color',cmap(trial-1,:));
    hold on;
    for i=2:40
        polarplot([0,prefL(i,trial)],[0,rNL(i,trial)],'color',cmap(trial-1,:))
        rlim([0,1])
        rticks([.5,1])
        set(gca,'fontsize',16)
    end
end
title('no-laser controls, trials 2-9')




subplot(2,4,5);hold on;
plot(prefL(:,1),prefL(:,2),'o','color',cmap(1,:))
plot([0,2*pi],[0,2*pi],'--k')
plot([0,2*pi],[pi,3*pi],'--k')
plot([0,2*pi],[-pi,pi],'--k')
xlim([0,2*pi]);
ylim([0,2*pi]);
set(gca,'fontsize',16)
xticks([0,pi,2*pi])
yticks([0,pi,2*pi])
xticklabels({'0','180','360'})
yticklabels({'0','180','360'})
xlabel('preferred arena heading (deg), trial 1')
ylabel('preferred arena heading (deg), trial 2')


subplot(2,4,6);hold on;
plot(prefNL(:,1),prefNL(:,2),'o','color',cmap(1,:))
plot([0,2*pi],[0,2*pi],'--k')
plot([0,2*pi],[pi,3*pi],'--k')
plot([0,2*pi],[-pi,pi],'--k')
xlim([0,2*pi]);
ylim([0,2*pi]);
set(gca,'fontsize',16)
xticks([0,pi,2*pi])
yticks([0,pi,2*pi])
xticklabels({'0','180','360'})
yticklabels({'0','180','360'})
xlabel('preferred arena heading (deg), trial 1')
ylabel('preferred arena heading (deg), trial 2')

subplot(2,4,7);hold on;
plot(prefNL(:,6),prefNL(:,7),'o','color',cmap(6,:))
plot([0,2*pi],[0,2*pi],'--k')
plot([0,2*pi],[pi,3*pi],'--k')
plot([0,2*pi],[-pi,pi],'--k')
xlim([0,2*pi]);
ylim([0,2*pi]);
set(gca,'fontsize',16)
xticks([0,pi,2*pi])
yticks([0,pi,2*pi])
xticklabels({'0','180','360'})
yticklabels({'0','180','360'})
xlabel('preferred arena heading (deg), trial 6')
ylabel('preferred arena heading (deg), trial 7')


subplot(2,4,8);hold on;
for i=1:8
    plot(prefNL(:,i),prefNL(:,i+1),'o','color',cmap(i,:))
end
plot([0,2*pi],[0,2*pi],'--k')
plot([0,2*pi],[pi,3*pi],'--k')
plot([0,2*pi],[-pi,pi],'--k')
xlim([0,2*pi]);
ylim([0,2*pi]);
set(gca,'fontsize',16)
xticks([0,pi,2*pi])
yticks([0,pi,2*pi])
xticklabels({'0','180','360'})
yticklabels({'0','180','360'})
xlabel('preferred heading (deg), trial i-1')
ylabel('preferred heading (deg), trial i')
end
