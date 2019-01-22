close all;
dt = 0.1;
rng = 0:dt:3;
sigma = 0.2;
gauss = dt/sqrt(2*pi())/sigma*exp(-0.5*((rng-mean(rng))/sigma).^2);
N = 30; % number of trials for comparison at the beginning and at the end of the session


for sess = 1;
    
    DDbef = [];
    DDaft = [];
    for k = 1:6
        D = gatherDrrd(k,sess,false);
        Dbef = D(1:N,1);
        Daft = D(end-N+1:end,1);
        DDbef = [DDbef; Dbef(:,1)];
        DDaft = [DDaft; Daft(:,1)];
    end
    DDbef1 = DDbef;
    DDaft1 = DDaft;
    
    
    DDbef = [];
    DDaft = [];
    for k = 7:15
        D = gatherDrrd(k,sess,false);
        Dbef = D(1:N,1);
        Daft = D(end-N+1:end,1);
        DDbef = [DDbef; Dbef(:,1)];
        DDaft = [DDaft; Daft(:,1)];
    end
    DDbef2 = DDbef;
    DDaft2 = DDaft;
    
    
    DDbef = [];
    DDaft = [];
    for k = [64 65 66 67 68 69]
        D = gatherDrrd(k,sess,false);
        Dbef = D(1:N,1);
        Daft = D(end-N+1:end,1);
        DDbef = [DDbef; Dbef(:,1)];
        DDaft = [DDaft; Daft(:,1)];
    end
    DDbef3 = DDbef;
    DDaft3 = DDaft;
    
    
    clf; hold on;
    
    %% plotting the beginning
    % plotting first group - first trials
    subplot(1,2,1);
    c = histc(DDbef1,rng-dt/2)/dt/length(DDbef1);
    g = conv(c,gauss,'same');
    plot(rng,g,'k-','linewidth',3);   

    % --- plotting the second group data
    c = histc(DDbef2,rng-dt/2)/dt/length(DDbef2);
    g = conv(c,gauss,'same');
    hold on;
    plot(rng,g,'-','color',[.6 .6 .6],'linewidth',3);
    
    % --- plotting the third group data
    c = histc(DDbef3,rng-dt/2)/dt/length(DDbef3);
    g = conv(c,gauss,'same');
    hold on;
    plot(rng,g,'-','color',[.8 .8 .8],'linewidth',3);
    
    % making up the figure
    ylim([0 1.4]);
    xlim([0 3]);
    set(gca,'box','on','fontsize',16)
    xlabel('t (s)');
    ylabel('probability density');
    title('Beginning');
    [h p] = kstest2(DDbef1,DDbef2);
    disp([h p]);
    
    
    %% plotting the end of the session
    subplot(1,2,2); hold on; 
    % plotting first group
    c = histc(DDaft1,rng-dt/2)/dt/length(DDaft1);
    g = conv(c,gauss,'same');
    plot(rng,g,'k-','linewidth',3);
    
    % second group
    c = histc(DDaft2,rng-dt/2)/dt/length(DDaft2);
    g = conv(c,gauss,'same');
    plot(rng,g,'-','color',[.6 .6 .6],'linewidth',3);
    
    % third group
    c = histc(DDaft3,rng-dt/2)/dt/length(DDaft3);
    g = conv(c,gauss,'same');
    plot(rng,g,'-','color',[.8 .8 .8],'linewidth',3);

    
    % making up figure 
    ylim([0 2.4]);
    set(gca,'yticklabel',{});
    ylim([0 1.4]);
    xlim([0 3]);
    xlabel('t (s)');
    legend({'timeout','no timeout' '1.2 s'});
    set(gca,'box','on','fontsize',16)
    title('End');
    [h p] = kstest2(DDaft1,DDaft2);
    disp([h p]);
    
end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
text(.015,0.95,'a','fontsize',20)
text(.500,0.95,'b','fontsize',20)


print('-depsc', 'fig3.eps');