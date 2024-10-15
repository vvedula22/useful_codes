close all; clear all; clc;

fname = "Figure11A\log_TONG.txt";
data = load(fname);
n  = size(data,2);
t  = data(:, 1);
V  = data(:, 2);
Ca_i = data(:, 3);
F = data(:, 4);

figure;
% subplot(311)
plot(t, V)
% ylabel("V (mV)")
xlabel("Time (ms)")
% subplot(312)
% plot(t, Ca_i)
% subplot(313)
% plot(t, (1000 < t & t < 3000) * -0.5)
