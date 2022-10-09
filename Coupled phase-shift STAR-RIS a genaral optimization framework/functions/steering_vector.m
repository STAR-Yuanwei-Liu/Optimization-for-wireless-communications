function [a] = steering_vector(array,theta, phi)
%STEERING_VECTOR Summary of this function goes here
%   Detailed explanation goes here
a = exp(-1i*array*K(theta,phi));
end

function k = K(theta,phi)
k = pi * [cos(theta).*cos(phi), sin(theta).*cos(phi), sin(phi)]';
end

