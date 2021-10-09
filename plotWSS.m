function plotWSS()
    close all;
    clear; clc;

    timeStepSize  = 1.0e-2;
    nTimePerCycle = 8000;
    wallFace = 'cyl_wall';
    
    [WSSdata, faceName, faceArea] = ...
        loadFluxData('results_test_strong/B_NS_WSS_average.txt');
    nFace  = size(WSSdata,2);
    nEnd   = size(WSSdata,1);
    %nStart = 1;
    nStart = nEnd -1*nTimePerCycle + 1;
    time = linspace(0, nEnd*timeStepSize, nEnd)';
    
    for iwall=1:nFace
        if strcmp(faceName(iwall),wallFace)
            break;
        end
    end
    if iwall > nFace
        error("Error: wall face not found")
    end
    WSS1 = WSSdata(:,iwall);
    
    [WSSdata, ~, ~] = ...
        loadFluxData('results_test_weak/B_NS_WSS_average.txt');
    WSS2 = WSSdata(:,iwall);

%     [Pdata, ~, ~] = ...
%         loadFluxData('results_test_strong/B_NS_Pressure_average.txt');
%     P1 = Pdata(:,iwall);
% 
%     [Pdata, ~, ~] = ...
%         loadFluxData('results_test_strong/B_NS_Pressure_average.txt');
%     P2 = Pdata(:,iwall);
     
    titleFont = 20;
    lineWidth = 2.5;
    axesFont  = 18;
    axesWidth = 1.5;
    
    %% Figure 1 - WSS vs. time
    figure('units','normalized','outerposition',[0.05 0.08 0.25 0.4]);
        plot(time(nStart:nEnd), WSS1(nStart:nEnd), ...
            'k-', 'LineWidth', lineWidth);
        hold on;
        plot(time(nStart:nEnd), WSS2(nStart:nEnd), ...
            'r-.', 'LineWidth', lineWidth/2);
        hold off;
        set(gca,'FontSize',axesFont,'LineWidth',axesWidth);
        xlabel('time','FontSize',axesFont);
        ylabel('WSS','FontSize',axesFont);
        ylim([0.0601 0.08]);
        ax = gca;
        ax.YTick = 0.06:0.004:0.08;
        ax.XTick = 20:20:100;
        title('Cylinder surface WSS','FontSize',titleFont);
        legend('Strong','Weak','FontSize',axesFont,'Location','northeast');
        legend('boxoff');

%     %% Figure 2 - Pressure vs. time
%     figure('units','normalized','outerposition',[0.05 0.08 0.25 0.4]);
%         plot(time(nStart:nEnd), P1(nStart:nEnd), ...
%             'k-', 'LineWidth', lineWidth);
%         hold on;
%         plot(time(nStart:nEnd), P2(nStart:nEnd), ...
%             'r-.', 'LineWidth', lineWidth/2);
%         hold off;
%         set(gca,'FontSize',axesFont,'LineWidth',axesWidth);
%         xlabel('time','FontSize',axesFont);
%         ylabel('Pressure','FontSize',axesFont);
%         title('Cylinder surface pressure','FontSize',titleFont);
%         legend('Strong','Weak','FontSize',axesFont,'Location','northeast');
%         legend('boxoff');

end

function [fData, faceL, areaL] = loadFluxData(fname)
    fid = fopen(fname,'r');

    % Read face names
    tline = fgetl(fid);
    strL  = split(tline);
    nFace = 0;
    for i=1:size(strL,1)
        if strL(i)~=""
            nFace = nFace + 1;
            faceL(nFace) = string(strL(i));
        end
    end

    % Read areas of each face
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
