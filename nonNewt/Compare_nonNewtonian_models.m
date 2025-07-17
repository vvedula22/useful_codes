clear all; clc; close all;

ns = 100;
sr = logspace(-6, 6, ns)';
mu = zeros(ns,5);
mu_x = zeros(ns,5);

% Carreau-Yasuda model
mu_i = [0.022, 0.035, 0.0345, 0.037]';
mu_o = [0.22, 1.6, 0.56, 0.592]';
lam = [0.11, 8.2, 1.902, 3.313]';
a = [0.644, 0.64, 1.25, 2.0]';
n = [0.392, 0.2128, 0.22, 0.3568]';

% Cassons model
mu_i(5) = 0.22803;
mu_i(6) = 0.1739;

mu_o(5) = 0.3953;
mu_o(6) = 0.6125;

strL = ["CY (a)"; "CY (b)"; "CY (c)"; "CY (d)"; "Cass (a)"; "Cass (b)"];
for i=1:ns
    for j=1:6
        if j<=4
            mu(i,j) = mu_i(j) + (mu_o(j)-mu_i(j))*...
                (1+(lam(j)*sr(i))^a(j))^((n(j)-1)/a(j));
            mu_x(i,j) = (mu_o(j)-mu_i(j))*(n(j)-1)*(lam(j)^a(j))*...
                (sr(i)^(a(j)-1))*(1+(lam(j)*sr(i))^a(j))^((n(j)-1)/a(j)-1);
        else
            if (sr(i) < 0.5)
                s = 0.5;
            else
                s = sr(i);
            end
            mu(i,j) = (mu_o(j)/sqrt(s) + mu_i(j))^2;
            mu_x(i,j) = 2*mu_o(j)*(s^(-1.5))*(mu_o(j)/sqrt(s) + mu_i(j));
        end
    end
end

fid = fopen('visc_nonNewt.dat','w');
fprintf(fid,'Variables=s, CYa, CYb, CYc, CYd, Cassa, Cassb\n');
for i=1:ns
    fprintf(fid,'%.9f', sr(i));
    for j=1:6
        fprintf(fid,'   %.9f', mu(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);


%% Make plot

titleFont = 14;
axesFont  = 14;
fontName  = "TimesNewRoman";

figure('units','normalized','outerposition',[0.05 0.3 0.35 0.45]);
    plot(sr, mu(:,1), 'r-', 'Linewidth', 2.0, ...
        'DisplayName', strL(1));
    hold on;
    plot(sr, mu(:,2), 'b-.', 'Linewidth', 2.0, ...
        'DisplayName', strL(2));
    plot(sr, mu(:,3), 'k--', 'Linewidth', 2.0, ...
        'DisplayName', strL(3));
    plot(sr, mu(:,4), 'r-v', 'Linewidth', 1.0, ...
        'MarkerSize', 4, 'DisplayName', strL(4));
    plot(sr, mu(:,5), 'g-s', 'Linewidth', 1.0, ...
        'MarkerSize', 4, 'DisplayName', strL(5));
    plot(sr, mu(:,6), 'm-o', 'Linewidth', 1.0, ...
        'MarkerSize', 4, 'DisplayName', strL(6));
    hold off;

    set(gca, 'FontName', fontName, 'FontSize',...
        axesFont, 'LineWidth', 1);
    xlim([1e-6 1e6]);
    ylim([0.01 5]);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    set( gca, 'Box', 'on', 'TickDir'     , 'out', ...
        'TickLength'  , [.01 .01], ...
        'YGrid'       , 'on', ...
        'XGrid'       , 'on', ...
        'XColor'      , [0 0 0 ], ...
        'YColor'      , [0 0 0 ], ...
        'LineWidth'   , 1.0);
    xlabel('Shear rate (s^-1)', 'FontSize', axesFont);
    ylabel('Viscosity (dyn-s/cm^2)', 'FontSize', axesFont);

    lgd = legend;
    lgd.FontSize = axesFont;
    lgd.Location = 'northeast';
    lgd.Orientation = 'vertical';
    lgd.EdgeColor = 'w';
    legend('boxon');


