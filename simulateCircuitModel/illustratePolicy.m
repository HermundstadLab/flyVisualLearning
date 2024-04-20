function illustratePolicy(p,pplot,plotNum)
% ILLUSTRATEPOLICY generates different figures that schematize the
% fixed-form policy:
%
% INPUTS
%   p: structure containing simulation parameters
%   pplot:  structure containing plotting parameters
%   plotNum key:
%       1: multiplicative operation 
%       2: phase shifts and readout of policy 
%       3: illustration of PFL activity
%       4: goal weight update
%       5: changing goal weights
%

nx = p.N;
thetaF = linspace(0,2*pi,nx+1);
theta  = thetaF(1:end-1);
headingLabels = 1:nx;


%illustrative goal heading
wBase = cos(theta-theta(3)).^2 + cos(theta-theta(17)).^3;
wBase = normalize(wBase);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,plotNum) 
    %initialize goal weights    
    w0 = wBase;
    wF = wBase;wF(nx+1) = wF(1);
    rW = sum(wBase.*exp(1i.*theta));
    goalW = mod(angle(rW),2*pi);


    figure;set(gcf,'Position',[200 200 1500 600],'color','w')
    
    %initialize headings
    np = 6;
    
    %vary goal and compass headings
    indShiftGoal = [0,0,0,0,0,0];
    indShiftHeading = [-4,0,4,0,0,0];
    gains = [1,1,1,1,.66,.33];
    indHeadings = 17;
    for k=1:numel(indShiftGoal)
        w = gains(k).*circshift(w0,indShiftGoal(k));
        rW = sum(w.*exp(1i.*theta));
        m0 = abs(rW);
        goalW = mod(angle(rW),2*pi);
        [~,indGoal] = min((goalW-theta).^2);
        
        subplot(3,np,k);
        
        cmap = flipud(gray(numel(indHeadings)+2));
        polarplot([0,theta(indHeadings+indShiftHeading(k))],[0,1],'color','k','linewidth',pplot.lw);hold on;
        polarplot([0,theta(indGoal)],[0,gains(k)],'color',pplot.cG,'linewidth',pplot.lw)
        set(gca,'fontsize',16)
        ax = gca;
        ax.RTick = [.5,1];
        ax.ThetaTick = [0,90,180,270];
        ax.RTickLabel = {};
        ax.ThetaTickLabel = {};

        subplot(3,np,np+k);hold on;
        cH = centerHeadingNeuronFB(theta,indHeadings+indShiftHeading(k));cH(nx+1) = cH(1);
        fill([thetaF,fliplr(thetaF),thetaF(1)],[zeros(size(thetaF)),fliplr(min([cH;[w,w(1)]])), 0],[.8,.8,.8],'LineStyle','none')
        plot(thetaF,cH,'linewidth',pplot.lw,'color','k');
        plot(theta(indHeadings),0,'o','markersize',12,'color','k')

        subplot(3,np,np+k);hold on;
        plot(thetaF,[w,w(1)],'linewidth',pplot.lw,'color',pplot.cG)
        plot(theta(indGoal),   0,'o','markersize',12,'color',pplot.cG)
        ylabel('heading activity')
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        set(gca,'fontsize',16)

        subplot(3,np,2*np+k);hold on;
        rO = zeros(1,nx);
        for i=1:nx
            if k<=4
                rO(i) = outputNeuron(theta,circIndex(i+nx/2,nx),w0,nx,'center');
            else
                rO(i) = outputNeuron(theta,circIndex(i+nx/2,nx),w,nx,'center');
            end
        end
        rO(nx+1) = rO(1);
        bar(thetaF,rO,'FaceColor','none','EdgeColor','k')
        plot([theta(indGoal),theta(indGoal)],[0,1],'-','linewidth',pplot.lw,'color',pplot.cG)
        ylabel('output activity')
        ylim([0,.8]);
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        set(gca,'fontsize',16)
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,plotNum) 
    %initialize goal weights
    w0 = wBase;
    w = w0;
    w = normalize(w);
    w0 = normalize(w0);
    rW = sum(w.*exp(1i.*theta));
    goalW = mod(angle(rW),2*pi);
    [~,indGoal] = min((goalW-theta).^2);

    
    figure;set(gcf,'Position',[200 200 1600 900],'color','w')
    
    %initialize heading
    readouts = {'right','center','left'};
    colors = [pplot.cS;pplot.cF;pplot.cS];
    shifts = [nx/4,nx/2,-nx/4];
    phases = {'+90','180','-90'};
    indHeading = 21;
    for i=1:numel(readouts)
         
        indPhaseShift = circIndex(indHeading+shifts(i),nx);
        subplot(3,3,i);
        polarplot([0,theta(indHeading)],[0,1],'color','k','linewidth',pplot.lw);
        hold on;
        polarplot([0,theta(indGoal)],[0,1],'color',pplot.cG,'linewidth',pplot.lw)
        polarplot([0,theta(indPhaseShift)],[0,1],'color',colors(i,:),'linewidth',pplot.lw)
        set(gca,'fontsize',16)
        ax = gca;
        ax.RTick = [.5,1];
        ax.ThetaTick = [0,90,180,270];
        ax.RTickLabel = {};
        ax.ThetaTickLabel = {};
        title(['phase shift = ',phases{i}]);

        cH  = centerHeadingNeuronFB(theta,indHeading);cH(nx+1)  = cH(1);
        cPS = centerHeadingNeuronFB(theta,indPhaseShift);cPS(nx+1) = cPS(1);
        wF  = w;wF(nx+1)= w(1); 
        subplot(3,3,3+i);hold on;
        fill([thetaF,fliplr(thetaF),thetaF(1)],[zeros(size(thetaF)),fliplr(min([cPS;wF])), 0],[.8,.8,.8],'LineStyle','none')
        plot(thetaF,cPS,'linewidth',pplot.lw,'color',colors(i,:));
        plot(thetaF,wF,'linewidth',pplot.lw,'color',pplot.cG)
        plot(theta(indHeading),0,'o','markersize',12,'color','k')
        plot(theta(indGoal),   0,'o','markersize',12,'color',pplot.cG)
        plot(theta(indPhaseShift), 0,'o','markersize',12,'color',colors(i,:))
        ylabel('heading activity')
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        set(gca,'fontsize',16)

        subplot(3,3,6+i);hold on;
        rO = zeros(1,nx);
        for j=1:nx
            rO(j) = outputNeuron(theta,j,w0,nx,readouts{i});
        end
        rO(nx+1) = rO(1);
        bar(thetaF,rO,'FaceColor','none','EdgeColor',colors(i,:))
        plot([theta(indGoal),theta(indGoal)],[0,1],'-','linewidth',pplot.lw,'color','k')
        ylabel('output activity')
        ylim([0.1,.8])
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        set(gca,'fontsize',16)
       
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(3,plotNum) 
    %initialize goal weights
    w0 = wBase;
    w0 = normalize(w0);

    figure;set(gcf,'Position',[200 200 1600 300],'color','w')
    
    %initialize heading
    readouts = {'right','center','left'};
    colors = [pplot.cS;pplot.cF;pplot.cS];
    shifts = [nx/4,nx/2,-nx/4];
    indHeading = 17;
    rO = zeros(3,nx);
    for i=1:numel(readouts)
        for j=1:nx
            rO(i,j) = outputNeuron(theta,j,w0,nx,readouts{i});
        end
        rO(i,nx+1) = rO(i,1);
    end
    subplot(1,3,1);hold on;
    bar(thetaF,rO(1,:)-rO(3,:),'FaceColor','none','EdgeColor',colors(1,:))
    ylabel('output activity')
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    title('R-L PFL3 activity')

    subplot(1,3,2);hold on;
    bar(thetaF,rO(2,:),'FaceColor','none','EdgeColor',colors(2,:))
    ylabel('output activity')
    ylim([0.1,.8])
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    title('PFL2 activity')

    subplot(1,3,3);hold on;
    bar(thetaF,rO(3,:)-rO(1,:),'FaceColor','none','EdgeColor',colors(3,:))
    ylabel('output activity')
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    title('L-R PFL3 activity')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(4,plotNum) 
    
    %location of heading
    indHeading  = 8;
    
    %define heading bump
    rC = .5*cos(theta-theta(indHeading)) + .5;
    
    %initialize goal weights to be updated
    w = .125*(-.2*cos(theta-theta(1))-.4*cos(theta-theta(5)).^2+1)+.4;

    dw = relu(rC-w).*nonlin(1-w) - relu(w-rC).*nonlin(w);
    wnew1 = w+p.alphaW.*dw*40;

    dw = -relu(rC-w).*nonlin(w) + relu(w-rC).*nonlin(1-w);
    wnew2 = w+p.alphaW.*dw*40;
    
    %rescale heading bump (for visualization)
    ymin = min([w,wnew1,wnew2]);
    ymax = max([w,wnew1,wnew2]);
    rC = rC*(ymax-ymin);
    rC = rC+ymin-min(rC);
    
    figure;hold on;set(gcf,'Position',[200 200 600 800],'color','w')
    subplot(2,1,1);hold on;
    plot(theta,rC,'-k','linewidth',2)
    plot(theta,w,'--','linewidth',2,'color',pplot.cG)
    plot(theta,wnew1,'-','linewidth',2,'color',pplot.cG)
    xlim([0,2*pi])
    set(gca,'fontsize',16)
    xlabel('heading')
    ylabel('activity')
    
    
    subplot(2,1,2);hold on;
    plot(theta,rC,'-k','linewidth',2)
    plot(theta,w,'--','linewidth',2,'color',pplot.cG)
    plot(theta,wnew2,'-','linewidth',2,'color',pplot.cG)
    xlim([0,2*pi])
    set(gca,'fontsize',16)
    xlabel('heading')
    ylabel('activity')
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(5,plotNum) 
    fmax = 1.5;
    
    %initialize goal weights
    w = wBase;
    rW = sum(w.*exp(1i.*theta));
    goalW = mod(angle(rW),2*pi);
    [~,indGoal] = min((goalW-theta).^2);
    strengthW = abs(rW)./12;
    
    wF = w;
    wF(nx+1) = wF(1);
    figure;set(gcf,'Position',[200 200 1600 900],'color','w')
    
    %initialize heading
    indHeading = 9;
    rC = centerHeadingNeuronFB(theta,indHeading);
    rCnorm = sum(rC.*exp(1i.*theta));
    strengthNorm = abs(rCnorm)./12;
    nsteps = 12;
    
    nB = 3;
    cmapG = buildColormap(pplot.cG,[1,1,1],nsteps+nB);
    cmapF = buildColormap(pplot.cF,[1,1,1],nsteps+nB);
    cmapS = buildColormap(pplot.cS,[1,1,1],nsteps+nB);
    
    subplot(4,4,3);
    polarplot([0,theta(indHeading)],[0,1],'color','k','linewidth',pplot.lw);
    hold on;
    polarplot([0,theta(indGoal)],[0,strengthW./strengthNorm],'-' ,'color',cmapG(nB,  :),'linewidth',pplot.lw)
    set(gca,'fontsize',16)
    ax = gca;
    ax.RTick = [.5,1];
    ax.ThetaTick = [0,90,180,270];
    ax.RTickLabel = {};
    ax.ThetaTickLabel = {};
    
    subplot(4,4,4);
    polarplot([0,theta(indHeading)],[0,1],'color','k','linewidth',pplot.lw);
    hold on;
    polarplot([0,pi/2],[0,strengthW./strengthNorm],'-' ,'color',cmapG(nB,  :),'linewidth',pplot.lw)
    set(gca,'fontsize',16)
    ax = gca;
    ax.RTick = [.5,1];
    ax.ThetaTick = [0,90,180,270];
    ax.RTickLabel = {};
    ax.ThetaTickLabel = {};
        
    subplot(4,4,7);hold on;
    rH = centerHeadingNeuronFB(theta,indHeading);rH(nx+1) = rH(1);
    plot(thetaF,rH,'linewidth',pplot.lw,'color','k');
    plot(thetaF,wF,'-', 'linewidth',pplot.lw,   'color',cmapG(nB,:))
    plot(theta(indHeading),0,'o','markersize',12,'color','k')
    plot(theta(indGoal),0,'o','markersize',12,'color',cmapG(nB,:))
    ylabel('heading activity')
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    
    subplot(4,4,8);hold on;
    wS = circshift(w,17-indGoal);
    rWS = sum(wS.*exp(1i.*theta));
    goalWS = mod(angle(rWS),2*pi);
    [~,indGoalS] = min((goalWS-theta).^2);
    
    wS(end+1) = wS(1);
    plot(thetaF,wS,'-', 'linewidth',pplot.lw,   'color',cmapG(nB,:))
    plot(theta(indGoalS),0,'o','markersize',12,'color',cmapG(nB,:))
    ylabel('heading activity')
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    
    nu = zeros(1,nx);
    for j=1:nx
        nu(j) = driftRate(theta,j,w,nx);
    end
    nuS = circshift(nu,17-indGoal);
    nuS(nx+1) = nuS(1);
    nu( nx+1) = nu(1);
    
    subplot(4,4,11);hold on;
    plot(thetaF,fmax./nu,'Color',cmapF(nB,:),'linewidth',pplot.lw)
    ylabel('output activity')
    xlim([0,2*pi])
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    
    subplot(4,4,12);hold on;
    plot(thetaF,fmax./nuS,'Color',cmapF(nB,:),'linewidth',pplot.lw)
    ylabel('output activity')
    xlim([0,2*pi])
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    
    
    pR = zeros(1,nx);
    for j=1:nx
        pR(j) = turnBias(theta,j,w,nx);
    end
    pRS = circshift(pR,17-indGoal);
    pRS(nx+1) = pRS(1);
    pR(nx+1)  = pR(1);
    
    subplot(4,4,15);hold on;
    plot(thetaF,pR,'Color',cmapS(nB,:),'linewidth',pplot.lw)
    ylabel('output activity')
    xlim([0,2*pi])
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
    
    subplot(4,4,16);hold on;
    plot(thetaF,pRS,'Color',cmapS(nB,:),'linewidth',pplot.lw)
    ylabel('output activity')
    xlim([0,2*pi])
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'fontsize',16)
        
        
    for i=1:nsteps    
        
        dw = relu(rC-w).*nonlin(1-w) - relu(w-rC).*nonlin(w); 
        wNew  = w+dw*.25;
        wNewF = wNew;
        wNewF(nx+1) = wNewF(1);
        rW = sum(wNew.*exp(1i.*theta));
        goalW = mod(angle(rW),2*pi);
        [~,indGoalNew] = min((goalW-theta).^2);
        strengthWnew = abs(rW)./12;
        
        subplot(4,4,3);
        polarplot([0,theta(indGoalNew)],[0,strengthWnew./strengthNorm],'-' ,'color',cmapG(i+nB,  :),'linewidth',pplot.lw)
        
        subplot(4,4,4);
        polarplot([0,pi/2],[0,strengthWnew./strengthNorm],'-' ,'color',cmapG(i+nB,  :),'linewidth',pplot.lw)
                
        if i<3
            
            subplot(4,4,i);
            polarplot([0,theta(indHeading)],[0,1],'color','k','linewidth',pplot.lw);
            hold on;
            polarplot([0,theta(indGoal)   ],[0,strengthW./strengthNorm],'--','color',cmapG(end-1,:),'linewidth',pplot.lw)
            polarplot([0,theta(indGoalNew)],[0,strengthWnew./strengthNorm],'-' ,'color',cmapG(end,  :),'linewidth',pplot.lw)
            set(gca,'fontsize',16)
            ax = gca;
            ax.RTick = [.5,1];
            ax.ThetaTick = [0,90,180,270];
            ax.RTickLabel = {};
            ax.ThetaTickLabel = {};
           

            subplot(4,4,4+i);hold on;
            rH = centerHeadingNeuronFB(theta,indHeading);rH(nx+1) = rH(1);
            plot(thetaF,rH,'linewidth',pplot.lw,'color','k');
            plot(thetaF,wF,  '--','linewidth',pplot.lw,'color',cmapG(end-1,:))
            plot(thetaF,wNewF,'-','linewidth',pplot.lw,'color',cmapG(end,:))
            plot(theta(indHeading),0,'o','markersize',12,'color','k')
            plot(theta(indGoal),   0,'o','markersize',12,'color',cmapG(end-1,:))
            plot(theta(indGoalNew),0,'o','markersize',12,'color',cmapG(end,:))
            ylabel('heading activity')
            xticks(0:pi/2:2*pi)
            xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
            set(gca,'fontsize',16)
            

        end
        
        subplot(4,4,7);hold on;
        plot(thetaF,wNewF,'-', 'linewidth',pplot.lw,'color',cmapG(i+nB,:))
        plot(theta(indGoalNew),0,'o','markersize',12,'color',cmapG(i+nB,:))
        ylabel('heading activity')
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        set(gca,'fontsize',16)
        
        subplot(4,4,8);hold on;
        wS = circshift(wNewF(1:end-1),17-indGoalNew);
        rWS = sum(wS.*exp(1i.*theta));
        goalWS = mod(angle(rWS),2*pi);
        [~,indGoalS] = min((goalWS-theta).^2);
    
        wS(end+1) = wS(1);
        plot(thetaF,wS,'-', 'linewidth',pplot.lw,'color',cmapG(i+nB,:))
        plot(theta(indGoalS),0,'o','markersize',12,'color',cmapG(i+nB,:))
        ylabel('heading activity')
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        set(gca,'fontsize',16)

    
        nu    = zeros(1,nx);
        nuOld = zeros(1,nx);
        for j=1:nx
            nu(j)    = driftRate(theta,j,wNew,nx);
            nuOld(j) = driftRate(theta,j,w,   nx);
        end
        nuS = circshift(nu,17-indGoalNew);
        nuS(nx+1)   = nuS(1);
        nu(nx+1)    = nu(1);
        nuOld(nx+1) = nuOld(1);
        
        
        if i<3
            subplot(4,4,8+i);hold on;
            plot(thetaF,fmax./nu,'Color',cmapF(end,:),'linewidth',pplot.lw)
            plot(thetaF,fmax./nuOld,'--','Color',cmapF(end-1,:),'linewidth',pplot.lw)
            ylabel('output activity')
            xticks(0:pi/2:2*pi)
            xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
            set(gca,'fontsize',16)
            ylim([0,15])
        end
        
        subplot(4,4,11);hold on;
        plot(thetaF,fmax./nu,'Color',cmapF(i+nB,:),'linewidth',pplot.lw)
        
        subplot(4,4,12);hold on;
        plot(thetaF,fmax./nuS,'Color',cmapF(i+nB,:),'linewidth',pplot.lw)
        
        

        pR = zeros(1,nx);
        pRold = zeros(1,nx);
        for j=1:nx
            pR(j)    = turnBias(theta,j,wNew,nx);
            pRold(j) = turnBias(theta,j,w,nx);
        end
        pRS = circshift(pR,17-indGoalNew);
        pR(nx+1)    = pR(1);
        pRold(nx+1) = pRold(1);
        pRS(nx+1)   = pRS(1);

        if i<3
            subplot(4,4,12+i);hold on;
            plot(thetaF,pR,'Color',cmapS(end,:),'linewidth',pplot.lw)
            plot(thetaF,pRold,'--','Color',cmapS(end-1,:),'linewidth',pplot.lw)
            ylabel('output activity')
            ylim([.15,.85])
            xticks(0:pi/2:2*pi)
            xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
            set(gca,'fontsize',16)
        end
        
       
        subplot(4,4,15);
        plot(thetaF,pR,'Color',cmapS(i+nB,:),'linewidth',pplot.lw)
        ylim([.15,.85])
        
        subplot(4,4,16);
        plot(thetaF,pRS,'Color',cmapS(i+nB,:),'linewidth',pplot.lw)
        ylim([.15,.85])
        
        wF = wNewF;
        w = wNew;
        indGoal = indGoalNew;
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end









