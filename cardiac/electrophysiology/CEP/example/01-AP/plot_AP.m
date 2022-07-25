close all; clear all; clc;
fname = 'log_AP.txt';
data = load(fname);
m = size(data,2);
t = data(:,1);
V = data(:,2);
ep = data(:,m-1);
Ta = data(:,m);

figure;
    subplot(2,2,1); plot(t, V);
    subplot(2,2,2); plot(t, Ta);
    subplot(2,2,3); plot(t, ep);
    subplot(2,2,4); plot(V, ep,'o');
    
    