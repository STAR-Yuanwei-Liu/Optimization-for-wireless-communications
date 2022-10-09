function [rate, W, theta_t, theta_r, phas_diff_all, sum_rate_all] = algorithm_PDD(para, G, h)
%The PDD-based algorithm for solve the throughput maximization problem subject to the coupled phase-shift constraints.
%  [rate, P, theta_t, theta_r, phas_diff_all, sum_rate_all] = algorithm_PDD(para, G, h)
%Inputs:
%   para: structure of the initial parameters
%   G: BS-STAR channel
%   h: STAR-user channel
%Outputs:
%   W: obtained beamformers
%   theta_t: obtained transmission coefficients of STAR-RIS
%   theta_r: obtained reflection coefficients of STAR-RIS
%   phas_diff_all: the absolute phase shift difference at each iteration of PDD-based algorithm
%   sum_rate_all: the throughput at each iteration of PDD-based algorithm
%Date: 01/10/2022
%Author: Zhaolin Wang


%% random initialization
theta_t = sqrt(0.5)*exp(1i * 2*pi*rand(para.N,1));
theta_r = sqrt(0.5)*exp(1i * 2*pi*rand(para.N,1));
theta_rr = theta_r; theta_tt = theta_t;
W = randn(para.M, para.K);
W = W ./ sqrt(trace(W*W')) * sqrt(para.Pt);

% dual variables
lambda_t = zeros(para.N, 1);
lambda_r = zeros(para.N, 1);

%% algorithm paramater
rho = 1; % penalty term
c = 0.1; % reduction factor of rho
epsilon = 1e-4; % convergence criteria


%% PDD-based algorthm
phas_diff_all = [];
sum_rate_all = [];

[gamma, ~] = SINR(para, W, theta_t, theta_r, G, h);
sum_rate = sum(log2(1 + gamma));
sum_rate_all = [sum_rate_all, sum_rate];
phas_diff_all = [phas_diff_all abs(angle(theta_r) - angle(theta_t))];

disp('%%%%%%%%%%%%%%%%%%%%%%%%%% Outer Loop - 1 %%%%%%%%%%%%%%%%%%%%%%%%%%');
eta = 10;
for i = 1:200

    % optimize AL problem
    [W, theta_t, theta_r, theta_tt, theta_rr, sum_rate, phas_diff] = algorithm_BCD(para, W, theta_t, theta_r, theta_tt, theta_rr, lambda_t, lambda_r, rho, G, h);
    
    % calculate constraint violation
    [v] = constraint_violation(theta_t, theta_r, theta_tt, theta_rr); 
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%% Outer Loop - ' num2str(i+1) ', Violation - ' num2str(v) ' %%%%%%%%%%%%%%%%%%%%%%%%%%']);
    sum_rate_all = [sum_rate_all, sum_rate];
    phas_diff_all = [phas_diff_all, phas_diff];
    
    if v < epsilon % algorithm converged
        break; 
    end

    if v <= eta    
        % update dual variables
        lambda_r = lambda_r + 1/rho*(theta_rr - theta_r);
        lambda_t = lambda_t + 1/rho*(theta_tt - theta_t);
    else
        % update penalty term
        rho = c*rho; 
        disp(['Penalty factor - ' num2str(rho)]);
    end

    eta = 0.9*v;
end

[gamma, ~] = SINR(para, W, theta_t, theta_r, G, h);
rate = log2(1 + gamma);
end


function [h] = constraint_violation(theta_t, theta_r, theta_tt, theta_rr)
    h1 = abs(theta_rr - theta_r); h1 = max(h1);
    h2 = abs(theta_tt - theta_t); h2 = max(h2);
    h = max([h1,h2]);
end

