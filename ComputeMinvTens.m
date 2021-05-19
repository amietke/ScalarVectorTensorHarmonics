function [MinvTens, Mb_tt, Mb_tp] = ComputeMinvTens(PSI_tt_lm,PHI_tt_lm,sin_t)
% Compute a generalized pseudo-inverse for least square harmonic projections
% of nematic tensor fields on a sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Total number of modes used
    N_modes = size(PSI_tt_lm,2);
    
    % Total number of grid points
    n_grid = size(PSI_tt_lm,1);
    
    % Sub-matrices
    M_Psi = 2 * ( PSI_tt_lm(2:end-1,:)'*PSI_tt_lm(2:end-1,:) ...
                + PHI_tt_lm(2:end-1,:)'*PHI_tt_lm(2:end-1,:) );
    M_Phi = M_Psi;
    M_PsiPhi = 2 * ( PSI_tt_lm(2:end-1,:)'*PHI_tt_lm(2:end-1,:) ...
                   - PHI_tt_lm(2:end-1,:)'*PSI_tt_lm(2:end-1,:));
    
    % % Build complete matrix, eliminate modes l = 0,1 for which 
    % no nematic tensor harmonic exists           
    M = [ M_Psi(5:end,5:end), M_PsiPhi(5:end,5:end); ...
          M_PsiPhi(5:end,5:end)', M_Phi(5:end,5:end) ];
    MinvTens = M^(-1);
    
    % Operators to turn numerical values of tensor Q_tt, Q_tp into correct
    % source vector to perform least-square approximation excluding poles
    Mb_tt = 2 * [ PSI_tt_lm(2:end-1,5:end)', zeros([N_modes-4, n_grid-2]); ...
                  zeros([N_modes-4, n_grid-2]), PHI_tt_lm(2:end-1,5:end)' ];
    
    PHI_tt_sint_lm = PHI_tt_lm(2:end-1,:).*repmat(1./sin_t,[1,N_modes]);
    PSI_tt_sint_lm = PSI_tt_lm(2:end-1,:).*repmat(1./sin_t,[1,N_modes]);
    
    Mb_tp = 2 * [ -PHI_tt_sint_lm(:,5:end)', zeros([N_modes-4, n_grid-2]); ...
                   zeros([N_modes-4, n_grid-2]), PSI_tt_sint_lm(:,5:end)' ];    
end

