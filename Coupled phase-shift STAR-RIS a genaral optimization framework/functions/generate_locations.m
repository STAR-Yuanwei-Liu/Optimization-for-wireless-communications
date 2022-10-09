function [BS_loc,user_loc,target_loc] = generate_locations(para)
%GENERATE_LOCATIONS Summary of this function goes here
%   Detailed explanation goes here

% in the format [distance, azimuth angle, elevation angle]

STAR_loc = [0,0,0];
BS_loc = [100, 40, 0];

% user locations
range = [50, 100];
user_loc = zeros(para.K, 3);
for i = 1:para.K
    user_loc(i,:) = [ (range(2)-range(1))*rand(1) + range(1), -180*rand(1), 180*rand(1)-90 ];
end

% target locations
target_loc = [30, 120, 30];

end

