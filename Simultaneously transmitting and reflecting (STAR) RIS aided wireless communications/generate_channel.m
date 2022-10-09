function [G, h] = generate_channel(para, BS_array, STAR_array)
%Generate the BS-user, RIS-user and BS-RIS channels 
%  [G, h] = generate_channel(para, BS_loc, user_loc, BS_array, STAR_array)
%Inputs:
%   para: structure of the initial parameters
%   BS_array: antenna array delopyment of BS
%   STAR_array: antenna array delopyment of STAR-RIS
%Outputs:
%   G: BS-RIS channel
%   h: STAR-user channels
%Date: 01/10/2022
%Author: Zhaolin Wang

epsilon = para.rician; % Rician factor

BS_loc = para.BS_loc;
user_loc = para.user_loc;

%% BS to STAR-RIS channel
% NLOS
G_NLOS = 1/sqrt(2) .* ( randn(para.N,para.M) + 1i*randn(para.N,para.M) );

% LOS
a_BR = steering_vector(BS_array, -BS_loc(2), -BS_loc(3));
a_RB = steering_vector(STAR_array, BS_loc(2), BS_loc(3));
G_LOS = a_RB*a_BR.';

% pathloss
path_loss = para.pathloss(BS_loc(1))';
path_loss = sqrt(10.^(- path_loss/10));
G = path_loss .* (sqrt(epsilon/(epsilon+1)) * G_LOS + sqrt(1/(epsilon+1)) * G_NLOS);


%% STAR-RIS to users channel
% NLOS
h_NLOS = 1/sqrt(2) .* ( randn(para.N,para.K) + 1i*randn(para.N,para.K) );

% LOS
h_LOS = zeros(para.N,para.K);
for k = 1:para.K
    h_LOS(:,k) = steering_vector(STAR_array, user_loc(k,2), user_loc(k,3));    
end

% pathloss
path_loss = para.pathloss(user_loc(:,1))';
path_loss = sqrt(10.^( (-para.noise_dB-path_loss)/10));
h = path_loss .* (sqrt(epsilon/(epsilon+1)) * h_LOS + sqrt(1/(epsilon+1)) * h_NLOS);


end