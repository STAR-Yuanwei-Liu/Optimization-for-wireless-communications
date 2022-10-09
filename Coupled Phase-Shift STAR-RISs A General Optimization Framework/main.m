clc
clear all
close all

addpath("functions\");
%% Parameters initialization
cvx_solver mosek
para = para_init();

%% Generate channels
[BS_array, STAR_array] = generate_arrays(para);
[G, h] = generate_channel(para, BS_array, STAR_array);

%% PDD-based algorithm
[rate, P, theta_t, theta_r, phas_diff_all, sum_rate_all] = algorithm_PDD(para, G, h);

%% Convergence
figure; 
subplot(2,1,1);
t = length(sum_rate_all);
plot(1:t, sum_rate_all, '-bo');
xlabel('Number of cumulative BCD iterations');
ylabel('Throughput (bit/s/Hz)');

subplot(2,1,2);
plot(phas_diff_all');
xlabel('Number of cumulative BCD iterations');
ylabel('Phase shift difference');
ylim([0,2*pi]);
yticks([0, pi/2 3*pi/2 2*pi]);
yticklabels({'0','\pi/2','3\pi/2' '2\pi'});