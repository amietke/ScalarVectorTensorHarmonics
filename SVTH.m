function [Y_lm_scal, Y_lm_vec, PSI_lm_vec, PHI_lm_vec, e_r, e_t, e_p, PSI_t_lm, PSI_p_lm, PHI_t_lm, PHI_p_lm, ...
            PSI_tt_lm, PSI_pp_lm, PSI_pt_lm, PHI_tt_lm, PHI_pp_lm, PHI_pt_lm] ...
                = SVTH(coord_grid,l_max)
% Computation of REAL-VALUED scalar, vector and tensor spherical 
% harmonics in on  the unit sphere
%
% INPUT:    coord_grid - Array with two columns of spherical coordinate angles
%                        [azimuthal angle phi, polar angle tgeta]. Grid points 
%                        can be  arbitrary, but first row must be [0,0] and 
%                        last row must be [0,pi] (the coordinate poles)
%           l_max - Compute all fields up to the l-th harmonic mode
%    
% OUTPUT:  Y_lm_scal - Scalar harmonics [(phi,theta), (l,m)]
%
%          Y_lm_vec, PSI_lm_vec, PHI_lm_vec - 
%          Cartesian representation of vector spherical harmonics in 
%          three-dimensional arrays [(phi,theta), (l,m), (x,y,z)]
%
%          e_r, e_t, e_p - Normalized spherical coordinate basis vectors
%                          [(phi,theta), (x,y,z)]
%
%          PSI_t_lm, PSI_p_lm, PHI_t_lm, PHI_p_lm - Fully covariant 
%          harmonic vector components [(phi,theta), (l,m)]
%
%          PSI_tt_lm, PSI_pp_lm, PSI_pt_lm, PHI_tt_lm, PHI_pp_lm, PHI_pt_lm -
%          Fully covariant component of the two trace-free symmetric 
%          spherical harmonic tensors [(phi,theta), (l,m)]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    N_modes = (l_max+1)^2;   % Total number of modes for given maximal l mode
    C_lm = zeros(N_modes,1); % Spherical harmonics normalizations
    
    lm_grid = length(coord_grid(:,1));
    phi = coord_grid(:,1);
    theta = coord_grid(:,2);
    cos_t = cos(theta);

    % Harmonic functions
    Y_lm_scal = zeros(lm_grid,N_modes);
    
    % Harmonic basis vectors cartesion components as VEC[(theta,phi), (l,m), (x,y,z)];
    Y_lm_vec = zeros(lm_grid,N_modes,3); 
    PSI_lm_vec = zeros(lm_grid,N_modes,3);
    PHI_lm_vec = zeros(lm_grid,N_modes,3);
    
    % Harmonic basis vectors spherical components VEC[(theta,phi), (l,m)]
    % (all lower indices)
    PSI_t_lm = zeros(lm_grid,N_modes);
    PSI_p_lm = zeros(lm_grid,N_modes);
    PHI_t_lm = zeros(lm_grid,N_modes);
    PHI_p_lm = zeros(lm_grid,N_modes);
    
    % Harmonic basis tensors TENS[(theta,phi), (l,m)] (all lower indices)
    PSI_pp_lm = zeros(lm_grid,N_modes);
    PSI_pt_lm = zeros(lm_grid,N_modes);
    PSI_tt_lm = zeros(lm_grid,N_modes);
    PHI_pp_lm = zeros(lm_grid,N_modes);
    PHI_pt_lm = zeros(lm_grid,N_modes);
    PHI_tt_lm = zeros(lm_grid,N_modes);
    
    % Unit vectors spherical coordinates
    e_r = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
    e_t = [cos(theta).*cos(phi), cos(theta).*sin(phi), -sin(theta)];
    e_p = [-sin(phi), cos(phi), zeros(lm_grid,1)];
    
    for l = 0:l_max % loop over l-modes
        Leg_l = legendre(l,cos_t); % P_l^m(x)
        if l == 0
            dLeg_l = zeros(1,lm_grid);
        else
            dLeg_l = legendre_derivative(l,cos_t); % dP_l^m(x)/dx
        end
        
        %%%% Prepare derivatives of Legendre polynomials as needed later %%%%
        dLeg_lm = zeros(lm_grid,l+1);  % dP_l^m[cos(t)]/dt
        cot_P_lm = zeros(lm_grid,l+1); % cot(t)*P_l^m[cos(t)] 
        Plm_sint = zeros(lm_grid,l+1); % P_l^m[cos(t)]/sin(t)
        for m = 0:l            
            % dP_l^m[cos(t)]/dt, cot(t)*P_l^m and P_l^m/sin(t) at poles
            dP_lm_0 = 0; % dP_l^m[cos(t)]/dt at theta = 0 (x=1)
            dP_lm_pi = 0; % dP_l^m[cos(t)]/dt at theta = pi (x=-1)
            cot_P_lm_0 = 0; % cot(t)*P_l^m[cos(t)] at theta = 0 (x=1)
            cot_P_lm_pi = 0; % cot(t)*P_l^m[cos(t)] at theta = pi (x=-1)
            if (m == 0)
                if l>0
                    % dP_l^(-1)[cos(t)]/dt is needed
                    P_l_leg_m1 = -Leg_l(2,:)/(l*(l+1));
                    dP_lm_0 = l*(l+1)*P_l_leg_m1(1) - Leg_l(2,1); %at theta = 0
                    dP_lm_pi = l*(l+1)*P_l_leg_m1(end) - Leg_l(2,end); % theta = pi
                end
            else
                for j = (abs(m)-1):min(abs(m)+1,l) 
                    if j == (abs(m)-1)
                        dP_lm_0 = (l+abs(m))*(l-abs(m)+1)*Leg_l(j+1,1);  % theta = 0
                        dP_lm_pi = (l+abs(m))*(l-abs(m)+1)*Leg_l(j+1,end);  % theta = pi
                        cot_P_lm_0 = dP_lm_0; % theta = 0
                        cot_P_lm_pi = dP_lm_pi; % theta = pi
                    end
                    if j == (abs(m)+1)
                        dP_lm_0 = dP_lm_0 - Leg_l(j+1,1); % theta = 0
                        dP_lm_pi = dP_lm_pi - Leg_l(j+1,end); % theta = pi
                        cot_P_lm_0 = cot_P_lm_0 + Leg_l(j+1,1); % theta = 0
                        cot_P_lm_pi = cot_P_lm_pi + Leg_l(j+1,end); % theta = pi
                    end
                end
                
            end
            
            % Full dP_l^m[cos(t)]/dt with correct pole values
            dP_lm_0 = -dP_lm_0/2;
            dP_lm_pi = -dP_lm_pi/2;
            dLeg_lm(:,m+1) = [dP_lm_0; -sin(theta(2:end-1)).*dLeg_l(abs(m)+1,2:end-1)'; dP_lm_pi];
            
            % Full cot(t)*P_l^m[cos(t)] and P_l^m[cos(t)]/sin(t) with finite pole values
            if m ~= 0                
                cot_P_lm_0 = -cot_P_lm_0/(2*abs(m));
                cot_P_lm_pi = -cot_P_lm_pi/(2*abs(m));
                cot_P_lm(:,m+1) = [cot_P_lm_0; cot(theta(2:end-1)).*Leg_l(abs(m)+1,2:end-1)'; cot_P_lm_pi];                
                Plm_sint(:,m+1) = [cot_P_lm_0; Leg_l(abs(m)+1,2:end-1)'./sin(theta(2:end-1)); -cot_P_lm_pi];
            end           
        end
        
        %%%% Calculate d^2 P_l^m(cos(t))/dt^2 and sin(t)^(-1)*(d/dt-cot(t))P_l^m %%%%
        ddLeg_lm = zeros(lm_grid,l+1); % d^2 P_l^m(cos(t))/dt^2
        cot_dLeg_lm = zeros(lm_grid,l+1); % cot(theta)*dP_l^m(cos(t))/dt
        diagel_Leg_lm = zeros(lm_grid,l+1); % sin(t)^(-1)*(d/dt-cot(t))P_l^m
        
        if l>0
            cot_P_l_m1 = -cot_P_lm(:,2)/(l*(l+1)); %cot(t)*P^(-1)_l(cos(t))
            
            % For calculation of sin(t)^(-1)*(d/dt-cot(t))P_l^m the
            % associated Legendre Polynomials of l-1 are needed
            Leg_l_lower = legendre(l-1,cos_t);
            
            for m = 0:l % New loop because it requires non-trivial terms across different m
                %%%%%%%%%%% Calculation of d^2 P_l^m(cos(t))/dt^2 %%%%%%%%%%%%%%%%
                
                % Prepare array cot(theta)*dP_l^m(cos(t))/dt   
                for j = (abs(m)-1):min(abs(m)+1,l) 
                    if j == (abs(m)-1)
                        if (m == 0)
                            cot_dLeg_lm(:,1) = -l*(l+1)*cot_P_l_m1/2;
                        else
                            cot_dLeg_lm(:,m+1) = -(l+abs(m))*(l-abs(m)+1)*cot_P_lm(:,j+1)/2;
                        end
                    end
                    if j == (abs(m)+1)
                        cot_dLeg_lm(:,m+1) = cot_dLeg_lm(:,m+1) + cot_P_lm(:,j+1)/2;
                    end
                end

                % Calculate d^2 P_l^m(cos(t))/dt^2 from recursion relation
                ddLeg_lm(:,m+1) = -(1-m)*cot_dLeg_lm(:,m+1) - m*Leg_l(m+1,:)';
                if m<l
                    ddLeg_lm(:,m+1) = ddLeg_lm(:,m+1) + cot_P_lm(:,m+2) + dLeg_lm(:,m+2);
                end                

                %%%%%% Calculation of sin(t)^(-1)*(d/dt-cot(t))P_l^m %%%%%%%%%%%            
                if m > 0 % Expression is only relevant for m ~= 0
                    if l > 1
                        P_llow_leg_m1 = -Leg_l_lower(2,:)'/((l-1)*l); % P^(-1)_(l-1)(cos(t))  
                    end                
                    for j = max(-1,1-l):min(m+2,l-1)
                        if j == m - 2
                            if j < 0
                                diagel_Leg_lm(:,m+1) = diagel_Leg_lm(:,m+1) ...
                                    + (l+m-2)*(l+m-1)*(l+m)*(l-m+1)*P_llow_leg_m1;
                            else
                                diagel_Leg_lm(:,m+1) = diagel_Leg_lm(:,m+1) ...
                                    + (l+m-2)*(l+m-1)*(l+m)*(l-m+1)*Leg_l_lower(j+1,:)';
                            end
                        end 
                        if j == m
                            diagel_Leg_lm(:,m+1) = diagel_Leg_lm(:,m+1) ...
                                            - 2*m*(l+m)*Leg_l_lower(j+1,:)';
                        end                        
                        if j == m + 2
                            diagel_Leg_lm(:,m+1) = diagel_Leg_lm(:,m+1) - Leg_l_lower(j+1,:)';
                        end
                    end
                    diagel_Leg_lm(:,m+1) = diagel_Leg_lm(:,m+1)/(4*m);
                end
            end
        end       

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    
        for m = -l:l % loop over all m-modes
            ind = l*(l+1) + m + 1; % Linear index through full mode-list            
            if m > 0
                N = sqrt(2)*cos(m*phi);
                dN = -sqrt(2)*m*sin(m*phi);
                ddN = -sqrt(2)*m^2*cos(m*phi);
            elseif m == 0
                N = ones(length(phi),1);
                dN = 0;
                ddN = 0;
            elseif m < 0
                N = sqrt(2)*sin(abs(m)*phi);
                dN = sqrt(2)*abs(m)*cos(abs(m)*phi);
                ddN = -sqrt(2)*abs(m)^2*sin(abs(m)*phi);
            end                       
            
            C_lm(ind) = (-1)^(1*abs(m)) * sqrt( (2*l+1)*factorial(l-abs(m)) / (4*pi*factorial(l+abs(m))) );
            Ylm_self = C_lm(ind) * Leg_l(abs(m)+1,:)' .* N;
            dtheta_Ylm = C_lm(ind) * dLeg_lm(:,abs(m)+1) .* N;
            ddtheta_Ylm = C_lm(ind) * ddLeg_lm(:,abs(m)+1) .* N;
            dphi_Ylm = C_lm(ind) * Leg_l(abs(m)+1,:)' .* dN;
            if m~= 0
                cot_dphi_Ylm = C_lm(ind) * cot_P_lm(:,abs(m)+1) .* dN; % cot(t)*dY_lm/dphi
                dphi_Ylm_sint = C_lm(ind) * Plm_sint(:,abs(m)+1) .* dN; % sin(t)^(-1)*dY_lm/dphi
            end
            ddphi_Ylm = C_lm(ind) * Leg_l(abs(m)+1,:)' .* ddN;
            if m~= 0
                ddphi_Ylm_sint = C_lm(ind) * Plm_sint(:,abs(m)+1) .* ddN; % sin(t)]^(-1)*d^2Y_lm/dphi^2
            end
            dphi_dtheta_Ylm = C_lm(ind) * dLeg_lm(:,abs(m)+1) .* dN;   % d^2Y_lm/(dphi*dt)
            diagel_Ylm = C_lm(ind) *  diagel_Leg_lm(:,abs(m)+1) .* dN; % sin(t)^(-1)*(d/dt-cot(t))*dY_lm/dphi            
            
            %%%% Normalized scalar harmonic function Y_{lm} %%%%
            Y_lm_scal(:,ind) = Ylm_self;

            %%%% Cartesian components of Spherical Vector Harmonics %%%%%
            Y_lm_vec(:,ind,1) = e_r(:,1).*Ylm_self;
            Y_lm_vec(:,ind,2) = e_r(:,2).*Ylm_self;
            Y_lm_vec(:,ind,3) = e_r(:,3).*Ylm_self;
            
            PSI_lm_vec(:,ind,1) = e_t(:,1).*dtheta_Ylm; 
            if m~=0
                PSI_lm_vec(:,ind,1) = PSI_lm_vec(:,ind,1) + e_p(:,1).*dphi_Ylm_sint;    
            end
            PSI_lm_vec(:,ind,2) = e_t(:,2).*dtheta_Ylm;
            if m~=0
                PSI_lm_vec(:,ind,2) = PSI_lm_vec(:,ind,2) + e_p(:,2).*dphi_Ylm_sint;
            end
            PSI_lm_vec(:,ind,3) = e_t(:,3).*dtheta_Ylm;
            
            PHI_lm_vec(:,ind,1) = e_r(:,2).*PSI_lm_vec(:,ind,3) - e_r(:,3).*PSI_lm_vec(:,ind,2);
            PHI_lm_vec(:,ind,2) = e_r(:,3).*PSI_lm_vec(:,ind,1) - e_r(:,1).*PSI_lm_vec(:,ind,3);
            PHI_lm_vec(:,ind,3) = e_r(:,1).*PSI_lm_vec(:,ind,2) - e_r(:,2).*PSI_lm_vec(:,ind,1);
            
            %%%% Spherical coordinate components of spherical vector harmonics %%%%
            PSI_t_lm(:,ind) = dtheta_Ylm;
            PSI_p_lm(:,ind) = dphi_Ylm;
            
            if m == 0
                PHI_t_lm(:,ind) = 0;
            else
                PHI_t_lm(:,ind) = -dphi_Ylm_sint;
            end
            PHI_p_lm(:,ind) = sin(theta).*dtheta_Ylm;         
            
            %%%% Spherical coordinate components of spherical tensor harmonics %%%%
            PSI_pp_lm(:,ind) = ddphi_Ylm + sin(theta).*cos(theta).*dtheta_Ylm ...
                                    + 0.5*l*(l+1)*sin(theta).^2.*Ylm_self; % Psi_pp
            if m == 0
                PSI_pt_lm(:,ind) = zeros(lm_grid,1);
            else
                PSI_pt_lm(:,ind) = dphi_dtheta_Ylm - cot_dphi_Ylm; % Psi_pt
            end
            PSI_tt_lm(:,ind) = ddtheta_Ylm + 0.5*l*(l+1)*Ylm_self; % Psi_tt = -Psi_p^p = -Psi_pp/[sin(t)^2]                        
            
            PHI_pp_lm(:,ind) = sin(theta).*dphi_dtheta_Ylm - cos(theta).*dphi_Ylm; % Phi_pp
            PHI_pt_lm(:,ind) = - cos(theta).*dtheta_Ylm ...
                                    - 0.5*l*(l+1)*sin(theta).*Ylm_self; % Phi_pt                                
            if m~=0
                PHI_pt_lm(:,ind) = PHI_pt_lm(:,ind) - ddphi_Ylm_sint; % ... - sin(t)^(-1)*d^2Y_lm/dphi^2
            end
            PHI_tt_lm(:,ind) = -diagel_Ylm; % Phi_tt = -Phi_p^p = -Phi_pp/[sin(t)^2] = -Psi_(tp)/sin(t)            
        end
         
        % Spherical tensor harmonics only exist for l>1 (remove numerical error)
        PSI_tt_lm(:,1:4) = zeros(lm_grid,4);
        PSI_pp_lm(:,1:4) = zeros(lm_grid,4);
        PSI_pt_lm(:,1:4) = zeros(lm_grid,4);
        PHI_tt_lm(:,1:4) = zeros(lm_grid,4);
        PHI_pp_lm(:,1:4) = zeros(lm_grid,4);
        PHI_pt_lm(:,1:4) = zeros(lm_grid,4);
    end
end