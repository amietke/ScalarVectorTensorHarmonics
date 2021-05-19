function MinvVec = ComputeMinvVec(PSI_t_nm,PHI_t_nm)
% Compute a generalized pseudo-inverse for least square harmonic projections
% of vector fields on a sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Sub-Matrices
    M_Psi = PSI_t_nm'*PSI_t_nm + PHI_t_nm'*PHI_t_nm;
    M_Phi = M_Psi;
    M_PsiPhi = PSI_t_nm'*PHI_t_nm - PHI_t_nm'*PSI_t_nm;

    % Build complete matrix, eliminate l = 0 for which no vector harmonic exists            
    M = [M_Psi(2:end,2:end), M_PsiPhi(2:end,2:end); ...
         M_PsiPhi(2:end,2:end)', M_Phi(2:end,2:end)];

    % Inverse          
    MinvVec = M^(-1);    
end

