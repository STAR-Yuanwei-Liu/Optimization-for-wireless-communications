function [BS_array, STAR_array] = generate_arrays(para)
%Generate the antenna array delpoyment of BS and STAR-RIS
%  [BS_array, STAR_array] = generate_arrays(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   BS_array: antenna array delopyment of BS
%   STAR_array: antenna array delopyment of STAR-RIS
%Date: 01/10/2022
%Author: Zhaolin Wang

BS_array = [(1:para.M)'-(para.M+1)/2, zeros(para.M,1), zeros(para.M,1)];

STAR_array = zeros(para.N, 3);
for i = 1:para.STAR_size(1)
    for  j = 1:para.STAR_size(2)
        n = (i-1)*para.STAR_size(2) + j;
        STAR_array(n,:) = [ i, 0, j ];
    end
end
STAR_array = STAR_array - [(para.STAR_size(1)+1)/2, 0, 0];

end

