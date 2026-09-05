close all; clear all; clc;

% Courtemanche-Ramirez-Nattel (1998) human atrial myocyte model
% Log columns (log_CRN.txt):
%   1  t (ms)
%   2  V (mV)      3  Na_i    4  K_i     5  Ca_i    6  Ca_up   7  Ca_rel
%   8  I_Na        9  I_K1   10  I_to   11  I_Kur  12  I_Kr   13  I_Ks
%  14  I_CaL      15  I_pCa  16  I_NaK  17  I_NaCa 18  I_bNa  19  I_bCa
%  20  I_rel      21  I_tr   22  I_up   23  I_up_leak
%  24  Fn         25  I_stim 26  I_sac
% Currents (cols 8-19, 25) are in absolute pA; divide by Cm = 100 pF
% for pA/pF.  SR fluxes (cols 20-23) are in mM/ms.

fname = 'log_CRN.txt';
data = load(fname);
t      = data(:,1);
V      = data(:,2);
Ca_i   = data(:,5);
Ca_up  = data(:,6);
Ca_rel = data(:,7);

Cm     = 100.0;
I_Na   = data(:,8) /Cm;
I_K1   = data(:,9) /Cm;
I_to   = data(:,10)/Cm;
I_Kur  = data(:,11)/Cm;
I_Kr   = data(:,12)/Cm;
I_Ks   = data(:,13)/Cm;
I_CaL  = data(:,14)/Cm;
I_NaK  = data(:,16)/Cm;
I_NaCa = data(:,17)/Cm;

I_rel  = data(:,20);
I_tr   = data(:,21);
I_up   = data(:,22);
Fn     = data(:,24);

figure('units','normalized','outerposition',[0.05 0.05 0.5 0.9]);

subplot(4,1,1);
plot(t, V, 'k-', 'LineWidth', 2);
ylabel('V (mV)');
title('Courtemanche (CRN) human atrial AP');
xlim([0 600]);

subplot(4,1,2);
plot(t, Ca_i*1e3, 'b-', 'LineWidth', 2);
ylabel('[Ca]_i (\muM)');
xlim([0 600]);

subplot(4,1,3);
plot(t, I_CaL, 'g-', t, I_to, 'r-', t, I_Kur, 'b-', 'LineWidth', 1.5);
ylabel('I (pA/pF)');
legend('I_{CaL}', 'I_{to}', 'I_{Kur}');
xlim([0 600]);

subplot(4,1,4);
plot(t, I_Kr, 'g-', t, I_Ks, 'r-', t, I_K1, 'b-', t, I_NaCa, 'm-', ...
    'LineWidth', 1.5);
ylabel('I (pA/pF)');
xlabel('t (ms)');
legend('I_{Kr}', 'I_{Ks}', 'I_{K1}', 'I_{NaCa}');
xlim([0 600]);
