function [values] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [values] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 01/10/2022
%Author: Zhaolin Wang

values.noise_dB = -110; % noise power in dBm
values.noise = 10^(values.noise_dB/10); 


values.M = 8; % overall antennas
values.STAR_size = [10,2]; % reflecting elements at RIS
values.N = values.STAR_size(1)*values.STAR_size(2); % reflecting elements at RIS

values.Pt = 10^(20/10); % overall transmit power
values.n = 1; % equivalent noise power
values.K = 4; % user number

values.pathloss= @(d) 30 + 22*log10(d); % path loss with d in m

values.rician = 10^(3/10); % rician factor

values.STAR_loc = [0,0,0]; %Location of STAR-RIS
values.BS_loc = [50, 20, 0]; % Location of BS

% user locations
range = [3,3];
values.user_loc = zeros(values.K, 3);
for i = 1:values.K/2
    values.user_loc(i,:) = [ (range(2)-range(1))*rand(1) + range(1), -180*rand(1), 180*rand(1)-90 ];
end

for i = values.K/2+1:values.K
    values.user_loc(i,:) = [ (range(2)-range(1))*rand(1) + range(1), 180*rand(1), 180*rand(1)-90 ];
end

end

