    close all;
    clear; clc;

    dt    = 1.0e-4;
    pconv = 1334.16;
    
    fid = fopen('flow.txt','r');
    C = textscan(fid,'%f %f');
    fData = cell2mat(C);
    fclose(fid);
    
    t = C{:,1};
    Q = -C{:,2};
    
    Rp = 3224.2321;
    C  = 1.09e-6;
    Rd = 12896.92865;
    
    Tc = t(end);
    nTime = int32(Tc/dt);
    nCycle = 10;
    
    ti = linspace(0,Tc,nTime)';
    Qi = interp1(t, Q, ti);
    
    P = zeros(nTime*nCycle,1);
    t = zeros(nTime*nCycle,1);
    k = 1;
    for n=1:nCycle
        for i=1:nTime
           t(k+1) = t(k) + dt;
           if i~=1 && i~=nTime
                dQdt = (Qi(i+1)-Qi(i-1))/(ti(i+1)-ti(i-1));
           elseif i == 1
                dQdt = (Qi(2)-Qi(1))/(ti(2)-ti(1));
           elseif i == nTime
                dQdt = (Qi(nTime)-Qi(nTime-1))/(ti(nTime)-ti(nTime-1));
           end
           r = (1/C) * (-P(k)/Rd + Qi(i)*(1+Rp/Rd)) + Rp*dQdt;

           P(k+1) = P(k) + dt*r;
           k = k + 1;
        end
    end
    P = P./pconv;
    
    titleFont = 20;
    lineWidth = 2.5;
    axesFont  = 18;
    axesWidth = 1.5;
    
    %% Figure 1
    figure('units','normalized','outerposition',[0.05 0.08 0.25 0.4]);
        plot(ti, Qi, 'k-', 'LineWidth', lineWidth);
        set(gca,'FontSize',axesFont,'LineWidth',axesWidth);
        xlabel('[sec]','FontSize',axesFont);
        ylabel('[mL/sec]','FontSize',axesFont);

    is = (nCycle-1)*nTime + 1;
    ie = nCycle*nTime;
    figure('units','normalized','outerposition',[0.05 0.5 0.25 0.4]);
        plot(t(is:ie)-t(is), P(is:ie), 'k-', 'LineWidth', lineWidth);
        set(gca,'FontSize',axesFont,'LineWidth',axesWidth);
        xlabel('[sec]','FontSize',axesFont);
        ylabel('[mmHg]','FontSize',axesFont);
