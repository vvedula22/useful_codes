close all; clear all; clc;

% Load data
fname  = 'log_PFIB.txt';
Tc     = 1000.0;  % ms
nTimePerCycle = 5000;

data   = load(fname);
n      = size(data,2);

% Data from last two cycles is plotted in figures
nend   = size(data,1);
nstart = nend - 2*nTimePerCycle + 1;

dt     = Tc / nTimePerCycle;

% Copy data to local variables
t      = data(nstart:nend,1);
V      = data(nstart:nend,2);
Ca_i   = data(nstart:nend,5);
I_Na   = data(nstart:nend,9);
I_K1   = data(nstart:nend,10);
I_to   = data(nstart:nend,11);
I_Kr   = data(nstart:nend,12);
I_Ks   = data(nstart:nend,13);
I_CaL  = data(nstart:nend,14);
I_NaCa = data(nstart:nend,15);
I_NaK  = data(nstart:nend,16);
I_fNa  = data(nstart:nend,26);
I_fK   = data(nstart:nend,27);
I_f    = I_fNa + I_fK;

% Reset time axis to 0
t(:) = t(:) - t(1) + dt;

%% Plots

f_AP = figure('units','normalized','outerposition',[0.65 0.2 0.3 0.3]);
plot(t./1000, V, 'k-', 'LineWidth', 4);
xlabel('Time (s)','FontSize',16,'FontWeight','bold');
ylabel('Action Potential (mV)','FontSize',16,'FontWeight','bold');
ylim([-80, 50]);
yticks([-80:40:40]);
set(gca,'FontSize',16,'FontWeight','bold');
movegui('east')

ylabel_pos = -0.13;

f_I = figure('units','normalized','outerposition',[0.05 0.1 0.3 0.9]);
figure(f_I);
    subplot(10,1,1,'Parent',f_I);
    plot(t./1000, I_Na, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{Na}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([-200, 0]);
    yticks([-200:200:0]);
    
    subplot(10,1,2,'Parent',f_I);
    plot(t./1000, I_CaL, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{CaL}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([-10, 0]);
    yticks([-10:10:0]);

    subplot(10,1,3,'Parent',f_I);
    plot(t./1000, I_to, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{to}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([0, 4]);
    yticks([0:4:4]);

    subplot(10,1,4,'Parent',f_I);
    plot(t./1000, I_Kr, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{Kr}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([0, 0.5]);
    yticks([0:0.5:0.5]);

    subplot(10,1,5,'Parent',f_I);
    plot(t./1000, I_Ks, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{Ks}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([-0.01, 0.4]);
    yticks([0:0.4:0.4]);

    subplot(10,1,6,'Parent',f_I);
    plot(t./1000, I_K1, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{K1}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([0, 0.25]);
    yticks([0:0.25:0.25]);

    subplot(10,1,7,'Parent',f_I);
    plot(t./1000, I_f, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{f}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([-0.1, 0.1]);
    yticks([-0.1:0.2:0.1]);

    subplot(10,1,8,'Parent',f_I);
    plot(t./1000, I_NaK, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{NaK}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([0.2, 0.4]);
    yticks([0.2:0.2:0.4]);

    subplot(10,1,9,'Parent',f_I);
    plot(t./1000, I_NaCa, 'k-', 'LineWidth', 3);
    ylh = ylabel(strjust(sprintf('I_{NaCa}\n(pApF^{-1})'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box','off','XTick',[],'FontSize',14,'FontWeight','bold');
    ylim([-0.5, 0.5]);
    yticks([-0.5:1.0:0.5]);

    subplot(10,1,10,'Parent',f_I);
    plot(t./1000, Ca_i.*1000, 'k-', 'LineWidth', 3);
    xlabel('Time (s)','FontSize',16,'FontWeight','bold');
    ylh = ylabel(strjust(sprintf('Ca^{2+}_i\n(uM)'),'center'),...
        'FontSize',14,'FontWeight','bold');
    ylh.Position(1) = ylabel_pos;
    set(gca,'Box', 'off', 'FontSize',14,'FontWeight','bold');
    ylim([0, 0.9]);
    yticks([0:0.9:0.9]);
movegui('west');

    
