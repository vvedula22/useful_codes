function plotFlux()
    close all;
    clear; clc;

    timeStepSize  = 1.1e-5;
    nTimePerCycle = 10000;
    pconv = 1334.16;
    
    [Qdata, faceName, faceArea] = ...
        loadFluxData('32-procs/B_NS_Velocity_flux.txt');

    [Pdata, ~, ~] = ...
        loadFluxData('32-procs/B_NS_Pressure_average.txt');
    Pdata = Pdata / pconv;

    nFace  = size(Qdata,2);
    nEnd   = size(Qdata,1);
    %nStart = 1;
    nStart = nEnd -1*nTimePerCycle + 1;
    time = linspace(0, nEnd*timeStepSize, nEnd)';
    
    titleFont = 20;
    lineWidth = 2.5;
    axesFont  = 18;
    axesWidth = 1.5;
    
    %% Figure 1
    figure('units','normalized','outerposition',[0.05 0.08 0.5 0.9]);
        m = floor(sqrt(nFace));
        n = ceil(nFace/m);
        for i=1:nFace
            subplot(m,n,i)
            plot(time(nStart:nEnd), Qdata(nStart:nEnd,i), ...
                'k-', 'LineWidth', lineWidth);

            set(gca,'FontSize',axesFont,'LineWidth',axesWidth);
            xlabel('[sec]','FontSize',axesFont);
            ylabel('[mL/sec]','FontSize',axesFont);
            %title(replace(faceName(i),'_','-'),'FontSize',titleFont);
            legend(replace(faceName(i),'_','-'), ...
                'FontSize', axesFont, 'Location', 'northeast');
            legend('boxoff');
        end
        sgtitle('Velocity Flux vs. time','FontSize',titleFont);

    %% Figure 1
    figure('units','normalized','outerposition',[0.05 0.08 0.5 0.9]);
        m = floor(sqrt(nFace));
        n = ceil(nFace/m);
        for i=1:nFace
            subplot(m,n,i)
            plot(time(nStart:nEnd), Pdata(nStart:nEnd,i), ...
                'k-', 'LineWidth', lineWidth);

            set(gca,'FontSize',axesFont,'LineWidth',axesWidth);
            %title('Average Pressure vs. time','FontSize',titleFont);
            xlabel('[sec]','FontSize',axesFont);
            ylabel('[mmHg]','FontSize',axesFont);
            legend(replace(faceName(i),'_','-'), ...
                'FontSize', axesFont, 'Location', 'northeast');
            legend('boxoff');
        end
        sgtitle('Average Pressure vs. time','FontSize',titleFont);
end

function [fData, faceL, areaL] = loadFluxData(fname)
    fid = fopen(fname,'r');

    % Read face names
    tline = fgetl(fid);
    faceL = string(strsplit(tline));

    % Read areas of each face
    nFace = size(faceL,2);
    areaL = zeros(nFace,1);
    tline = fgetl(fid);
    str   = string(strsplit(tline));
    areaL = str2double(str);

    % Rewind and read all data
    frewind(fid);
    fmt = repmat('%f ', 1, nFace);
    C = textscan(fid,fmt,'HeaderLines',3);
    fData = cell2mat(C);
    fclose(fid);
end
