close all; clear all; clc;

fname = "log_TONG.txt";
data = load(fname);
n    = size(data,2);
t    = data(:, 1);
V    = data(:, 2);
Ca_i = data(:, 3);
F    = data(:, n);

%% Plots
figure('units','normalized','outerposition',[0.65 0.4 0.6 0.4]);
    %% Voltage, Force
    subplot(1,2,1);
    plot(t./1000.0, V, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t./1000.0, F, 'r--', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 20.0];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [-75 25];
    ylabel('Action Potential (mV)', ...
        'FontSize',18,'FontWeight','bold','Color','k');
    yticks([-75:25:25]);
    
    yyaxis right;
    ax.YLim = [0.51 3.5];
    ylabel('Force (uN)', ...
        'FontSize',18,'FontWeight','bold','Color',[0.8 0 0]);
    yticks([0:0.5:4]);

    dim = [0.3 0.75 0.25 0.1];
    str = {'Tong USMC model'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;

    %% Ca_i, Force
    subplot(1,2,2);
    plot(t./1000.0, Ca_i.*1e6, 'k-', 'LineWidth', 4);
    hold on;
    yyaxis right;
    plot(t./1000.0, F, 'r--', 'LineWidth', 4);
    
    set(gca, 'FontName', 'Times', 'FontSize', 18, 'LineWidth', 2);
    ax = gca;
    ax.XLim = [0.0 20.0];
    xlabel('Time (s)',...
        'FontSize',18,'FontWeight','bold','Color','k');
    
    yyaxis left;
    ax.YLim = [1 400];
    ylabel('Ca$_i$ (nM)', 'Interpreter','latex', ...
        'FontSize',24,'FontWeight','bold','Color','k');
    yticks([0:100:400]);
    
     yyaxis right;
     ax.YLim = [0.51 3.5];
     ylabel('Force (uN)', ...
         'FontSize',18,'FontWeight','bold','Color',[0.8 0 0]);
    yticks([0:0.5:4]);

    dim = [0.75 0.75 0.25 0.1];
    str = {'Tong USMC model'};
    annotation('textbox',dim,'String',str,...
        'FitBoxToText','on','FontSize',20, 'FontWeight','bold',...
        'BackgroundColor','w','EdgeColor','w');
    hold off;


