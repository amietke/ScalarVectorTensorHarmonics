function HarmonicProjectionExamples
% This file provides examples of how to use the least-square harmonic
% projections files LSQscalar.m, LSQvector.m, LSQtensor.m by performing
% transformations between real-space fields and harmonic coefficients.
% In each example, a random field is generated, and its harmonic representation
% is recovered and compared to the ground truth.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Clear figures
    close all;    
    
    %%%%% SELECT MAXIMAL MODE NUMBER AND GRID FOR PROJECTION EXAMPLES %%%%% 
    % Maximal mode number used for these examples
    l_max = 4;
    
    %%%%%% Examples of grids on unit sphere (comment out as needed) %%%%%%
    %%%--- (1) Regular grid ---%%%
    % Angular resolution of spherical coordinate grid
    % (Should divide 360)
    p_res = 20; % Azimuthal angle resolution in degree
    t_res = 5; % Polar angle resolution in degree 
    
    % Azimuthal angle phi varies in [0,360], Polar angle theta varies in [0,180]
    coord_grid = res2grid(p_res,t_res); % [phi,theta] 
    
%     %%%--- (2) Random grid ---%%%
%     Nrand = 800;
%     rng('shuffle');
%     coord_grid = [0,0; 2*pi*rand([Nrand-2,1]), ...
%                   acos(2*rand([Nrand-2,1])-1); 0, pi];
    %%%%%%%%%%%% End grid example s%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    % Total number of modes for given maximal l mode
    N_modes = (l_max+1)^2;     
            
    % Calculate components of harmonic vector and tensor fields on this grid
    [Y_lm_scal, ~, PSI_lm_vec, PHI_lm_vec, ~, ~, ~, PSI_t_lm, ~, PHI_t_lm, ~, ...
     PSI_tt_lm, ~, PSI_tp_lm, PHI_tt_lm, ~, PHI_tp_lm] = SVTH(coord_grid,l_max);    
    
    %%%%%%%%%%%%%%% EXAMPLE 1: SCALAR harmonic projection %%%%%%%%%%%%%
    % Generate a random sequence of harmonic mode coefficients in [-1,1]
    f_lm = 2 * ( rand([N_modes,1]) - 0.5 );
    
    % Compute the corresponding scalar function in real space
    f_test = Y_lm_scal*f_lm;
    
    %-- Starting point in most practical applications: Some scalar function 
    %-- f_test is given on discrete points on the sphere. 
    %-- Goal: Find a representation in terms of scalar spherical harmonics
    
    % Generate relevant inversion operator
    MinvScal = pinv(Y_lm_scal);
    
    % Perform harmonic projection to find mode decomposition    
    flm_proj = LSQscalar(f_test,MinvScal);
    
    % Error of scalar least-sq projection
    Err_scalar = mean(abs(f_lm-flm_proj));
    display(Err_scalar);
    %%%%%%%%%%%%%%%%%%%%%%%%% END EXAMPLE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%% EXAMPLE 2: VECTOR harmonic projection %%%%%%%%%%%%%
    % Generate a random sequence of harmonic mode coefficients in [-1,1]
    v1_lm = [0; 2 * ( rand([N_modes-1,1]) - 0.5 )];
    v2_lm = [0; 2 * ( rand([N_modes-1,1]) - 0.5 )];    
    
    % Cartesian components of the corrsponding vector field 
    v_test = [ PSI_lm_vec(:,:,1)*v1_lm + PHI_lm_vec(:,:,1)*v2_lm, ... % x-comp.
               PSI_lm_vec(:,:,2)*v1_lm + PHI_lm_vec(:,:,2)*v2_lm, ... % y-comp.
               PSI_lm_vec(:,:,3)*v1_lm + PHI_lm_vec(:,:,3)*v2_lm ];   % z-comp.
    
    %-- Starting point in most practical applications: Cartesian components of
    %-- some vector field v_test are given on discrete points on the sphere. 
    %-- Goal: Find a representation in terms of vector spherical harmonics
    
    % Scalar product of Cartesian vector field with 
    % Cartesian representation of vector spherical harmonics
    PsiDotVec = PSI_lm_vec(:,:,1)'*v_test(:,1) ...
              + PSI_lm_vec(:,:,2)'*v_test(:,2) ...
              + PSI_lm_vec(:,:,3)'*v_test(:,3);
    PhiDotVec = PHI_lm_vec(:,:,1)'*v_test(:,1) ...
              + PHI_lm_vec(:,:,2)'*v_test(:,2) ...
              + PHI_lm_vec(:,:,3)'*v_test(:,3);
           
    % Generate relevant inversion operator
    MinvVec = ComputeMinvVec(PSI_t_lm,PHI_t_lm);
    
    % Perform harmonic projection to find mode decomposition    
    [v1lm_proj, v2lm_proj] = LSQvector(PsiDotVec,PhiDotVec,MinvVec);
    
    % Error of vector harmonics least-sq projection
    Err_vector1 = mean(abs(v1_lm-v1lm_proj));
    display(Err_vector1);
    
    Err_vector2 = mean(abs(v2_lm-v2lm_proj));
    display(Err_vector2);
    %%%%%%%%%%%%%%%%%%%%%%%%% END EXAMPLE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%% EXAMPLE 3: Nematic tensor harmonic projection %%%%%%%%
    % Generate a random sequence of harmonic mode coefficients in [-1,1]
    q1_lm = [ zeros([4,1]) ; 2 * ( rand([N_modes-4,1]) - 0.5 )];
    q2_lm = [ zeros([4,1]) ; 2 * ( rand([N_modes-4,1]) - 0.5 )];    
    
    % Fully covariant components with respect to a surface parameterization
    % in terms of spherical coordinates t...theta, p...phi
    Q_tt = PSI_tt_lm*q1_lm + PHI_tt_lm*q2_lm; % Covariant theta,theta component
    Q_tp = PSI_tp_lm*q1_lm + PHI_tp_lm*q2_lm; % Covariant theta,phi component
    
    %-- Starting point for trace-free symmetric tensor harmonics projections: 
    %-- Fully covariant tensor components (theta,theta) and (theta,phi) are 
    %-- given on some grid on the sphere
    %-- Goal: Find a representation in terms of nematic tensor spherical harmonics
    
    % Generate relevant inversion operators
    [MinvTens, Mb_tt, Mb_tp] = ...
        ComputeMinvTens( PSI_tt_lm,PHI_tt_lm,sin(coord_grid(2:end-1,2)) );
    
    % Perform harmonic projection on tensor components to find mode decomposition    
    [q1lm_proj, q2lm_proj] = LSQtensor(Q_tt,Q_tp,MinvTens,Mb_tt,Mb_tp);
    
    % Error of tensor harmonics least-sq projection
    Err_tens1 = mean(abs(q1_lm-q1lm_proj));
    display(Err_tens1);
    
    Err_tens2 = mean(abs(q2_lm-q2lm_proj));
    display(Err_tens2);
    %%%%%%%%%%%%%%%%%%%%%%%%% END EXAMPLE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%    
end