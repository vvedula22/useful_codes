close all; clear all; clc;

fname = 'log_NYG.txt';
data = load(fname);
t      = data(:,1);
V      = data(:,2);

I_Na   = data(:,19);
I_CaL  = data(:,20)
I_t    = data(:,21);
I_sus  = data(:,22);

I_K1   = data(:,25);
I_NaK  = data(:,28);
I_NaCa = data(:,30);

Ca_i   = data(:,8);

O_C     = data(:,14);
O_TC    = data(:,15);
O_TMgC  = data(:,16);
O_Calse = data(:,18);

I_up    = data(:,31);
I_rel   = data(:,33);

Ca_up   = data(:,10);
Ca_rel  = data(:,11);

%% Plots
figure('units','normalized','outerposition',[0.3 0.3 0.3 0.5]);
    subplot(3,1,1);
    plot(t, V, 'k-', 'LineWidth', 3);   
    box off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 0.5];
    ax.XTickLabel = [];
    ax.YLim = [-80 40];
    ylabel('mV','FontSize',20,'FontWeight','bold','Color','k');
    yticks([-80:40:40]);
    
    subplot(3,1,2);
    plot(t, I_Na,'k-','LineWidth', 3);
    box off;
    hold on;
    plot(t, I_CaL,'k--','LineWidth',3);
    plot(t, I_sus,'k-.','LineWidth',3);
    plot(t, I_t,'k.','LineWidth',3);
    hold off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 0.5];
    ax.XTickLabel = [];
    ax.YLim = [-300 400];
    ylabel('pA','FontSize',20,'FontWeight','bold','Color','k');
    yticks([-300:200:400]);

    subplot(3,1,3);
    plot(t, I_NaCa,'k-','LineWidth', 3);
    box off;
    hold on;
    plot(t, I_NaK,'k--','LineWidth', 3);
    plot(t, I_K1,'k-.','LineWidth', 3);
    hold off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 0.5];
    xlabel('Time (s)','FontSize',18,'FontWeight','bold','Color','k');
    ax.YLim = [-40 60];
    ylabel('pA','FontSize',20,'FontWeight','bold','Color','k');
    yticks([-40:20:60]);

%%
figure('units','normalized','outerposition',[0.3 0.3 0.3 0.75]);
    subplot(5,1,1);
    plot(t, V, 'k-', 'LineWidth', 3);   
    box off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 0.5];
    ax.XTickLabel = [];
    ax.YLim = [-80 40];
    ylabel('mV','FontSize',20,'FontWeight','bold','Color','k');
    yticks([-80:40:40]);
    
    subplot(5,1,2);
    plot(t, Ca_i.*1000, 'k-', 'LineWidth', 3);   
    box off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0 0.5];
    ax.XTickLabel = [];
    ax.YLim = [-0.1 1.4];
    ylabel('$\mu$mol/L','Interpreter','Latex',...
        'FontSize',20,'FontWeight','bold','Color','k');
    yticks([0:0.4:1.2]);

    subplot(5,1,3);
    plot(t, O_C,'k-','LineWidth', 3);
    box off;
    hold on;
    plot(t, O_TC,'k--','LineWidth',3);
    plot(t, O_TMgC,'k-.','LineWidth',3);
    plot(t, O_Calse,'k.','LineWidth',3);
    hold off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0 0.5];
    ax.XTickLabel = [];
    ax.YLim = [-0.1 0.6];
    yticks([0:0.2:0.4]);

    subplot(5,1,4);
    plot(t, I_rel./1000,'k-','LineWidth', 3);
    box off;
    hold on;
    plot(t, I_up./1000,'k--','LineWidth', 3);
    hold off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 0.5];
    xlabel('Time (s)','FontSize',18,'FontWeight','bold','Color','k');
    ax.YLim = [0 12];
    ylabel('nA','FontSize',20,'FontWeight','bold','Color','k');
    yticks([0:4:12]);

    subplot(5,1,5);
    plot(t, Ca_up,'k-','LineWidth', 3);
    box off;
    hold on;
    plot(t, Ca_rel,'k--','LineWidth', 3);
    hold off;
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 0.5];
    xlabel('Time (s)','FontSize',18,'FontWeight','bold','Color','k');
    ax.YLim = [0 1];
    ylabel('mmol/L','FontSize',20,'FontWeight','bold','Color','k');
    yticks([0:0.2:1]);

