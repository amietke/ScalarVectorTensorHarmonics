function f_lm = LSQscalar(f,Pinv_Ylm)
% Harmonic transformation of a scalar function in the least-squares sense
%
% INPUT:    f - scalar function values on surface grid points
%           Pinv_Ylm - Pseudo-inverse of scalar harmonic function array
%    
% OUTPUT:   f_lm - Harmonic coefficients 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f_lm = Pinv_Ylm*f;
end