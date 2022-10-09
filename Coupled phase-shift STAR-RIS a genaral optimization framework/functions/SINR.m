function [gamma, upsilon] = SINR(para, W, theta_t, theta_r, G, g)
%Calculate the signal-to-interference-plus-noise ratio (SINR) at each communication user
%  [G, h] = generate_channel(para, BS_loc, user_loc, BS_array, STAR_array)
%Inputs:
%   para: structure of the initial parameters
%   W: Beamformers
%   theta_t: transmission coefficients of STAR-RIS
%   theta_r: reflection coefficients of STAR-RIS
%Outputs:
%   gamma: SINR
%   upsilon: auxiliary in WMMSE
%Date: 01/10/2022
%Author: Zhaolin Wang

Theta_t = diag(theta_t); ht = G' * Theta_t' * g; 
Theta_r = diag(theta_r); hr = G' * Theta_r' * g;  

gamma = zeros(para.K,1);
upsilon = zeros(para.K,1);
for k = 1:para.K
    if k <= para.K/2
        hk = ht(:,k); 
    else
        hk = hr(:,k);
    end   
    wk = W(:,k);
    upsilon(k) = hk'*wk / (real(hk'*(W*W')*hk) + 1);
    W_inter = W;
    W_inter(:,k) = [];
    gamma(k) = abs(hk'*wk)^2 / (real(hk'*(W_inter*W_inter')*hk) + 1);   
end

end

