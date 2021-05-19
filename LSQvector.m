function [v1_lm, v2_lm] = LSQvector(PsiDotVec,PhiDotVec,MinvVec)
% Harmonic transformation of a vector with given
% numerical components v_i in the least-squares sense
%
% INPUT:    M_inv - inverse least-square matrix (see notes)
%           PsiDotVec,PhiDotVec: Scalar products \mathbf{v}cdot\boldsymbol{V}_{lm}
%           for harmonic basis vectors Psi_lm and Phi_lm
%    
% OUTPUT:   Vector of mode coefficients v_nm = [v^(1)_lm, v^(2)_lm]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Total number of modes used 
    N_modes = length(PsiDotVec);

    % Calculate source vector (eliminate l = 0 elements of dot-arrays)
    b = [PsiDotVec(2:end); PhiDotVec(2:end)];
    
    % Solve linear regression M*v_lm = b
    v_lm = MinvVec*b;
    
    % Arrange mode coefficients and fill in zeros for modes with l=0 and l=1    
    v1_lm = [0; v_lm(1:N_modes-1)];  % v^(1)_lm
    v2_lm = [0; v_lm(N_modes:end)];  % v^(2)_lm
end