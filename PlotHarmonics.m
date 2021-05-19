function PlotHarmonics
% Plot scalar, vector and trace-free sysmmetric tensor harmonics for a given mode (l,m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%% CHOOSE MODE AND PLOT SETTINGS %%%%%%%%%%%%%%%%%%
    % Select mode number for plotting
    l = 2;  % Harmonic mode l
    m = -1;  % Angular mode number ( |m| <= l !!!)
        
    % Angular resolution of spherical coordinate grid for plotting
    % (Should divide 360)
    p_res = 9; % Azimuthal angle resolution in degree
    t_res = 9; % Polar angle resolution in degree               
    
    % Nematic axis size scaling for plotting
    nem_sc = 0.06;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Clear figures
    close all;
    
    % Linear index through full mode-list
    ind_mode = l*(l+1) + m + 1; 
    
    % Azimuthal angle phi varies in [0,360], Polar angle theta varies in [0,180]
    coord_grid = res2grid(p_res,t_res); % [phi,theta]

    % Calculate components of harmonic vector and tensor fields on the grid
    [Y_lm_scal, ~, PSI_lm_vec, PHI_lm_vec, e_r, e_t, e_p, ~, ~, ~, ~, ...
     PSI_tt_lm, ~, PSI_pt_lm, PHI_tt_lm, ~, PHI_pt_lm] = SVTH(coord_grid,l);  
                    
    %%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%--- EXAMPLE 1: Scalar harmonic plotted on a sphere ---%%%%%%%%%
    Y_lm_scal_sel = Y_lm_scal(:,ind_mode); % Extract relevant harmonic
    ax1 = MakeSphereAxis(p_res,t_res);     % Generate surface plot template
    ax1.CData = Grid2Mesh(Y_lm_scal_sel,p_res,t_res); % Plot
    
    % Add plot title
    title(['$Y_{',num2str(l),',',num2str(m),'}$'],...
           'Interpreter','latex','Fontsize',30);
    
    %%%--- EXAMPLE 2: Scalar harmonic plotted on the Phi-Theta domain ---%%
    [P, T] = meshgrid(pi*(0:p_res:360)/180,pi*(0:t_res:180)/180);
    figure;
    pcolor(P,T,Grid2Mesh(Y_lm_scal_sel,p_res,t_res));
    
    % Add axis labels and plot title
    xlabel('$\theta$','Interpreter','latex','Fontsize',30);
    ylabel('$\phi$','Interpreter','latex','Fontsize',30,'rot',0);
    title(['$Y_{',num2str(l),',',num2str(m),'}$'],...
           'Interpreter','latex','Fontsize',30);
    
    %%%%--- EXAMPLE 3: Vector harmonic: (Psi_lm)_i for chosen mode ---%%%%%
    Psi_Cart = [ PSI_lm_vec(:,ind_mode,1), ...
                 PSI_lm_vec(:,ind_mode,2), ...
                 PSI_lm_vec(:,ind_mode,3) ]; % Cartesian vector representation
    
    ax2 = MakeSphereAxis(p_res,t_res);     % Generate surface plot template
    ax2.CData = Grid2Mesh(vecnorm(Psi_Cart,2,2),p_res,t_res); % Plot vector norm
    
    % Add vector field directions
    hold on
    quiver3(e_r(:,1), e_r(:,2), e_r(:,3), Psi_Cart(:,1), Psi_Cart(:,2), Psi_Cart(:,3),...
            'Color','r','Linewidth',1.5,'AlignVertexCenters','on','AutoScaleFactor',1.2);         
    
    % Add plot title
    title(['$\Psi_{i}^{(',num2str(l),',',num2str(m),')}$'],...
           'Interpreter','latex','Fontsize',30);    
    
    %%%%--- EXAMPLE 4: Vector harmonic: (Phi_lm)_i for chosen mode ---%%%%%
    Phi_Cart = [ PHI_lm_vec(:,ind_mode,1), ...
                 PHI_lm_vec(:,ind_mode,2), ...
                 PHI_lm_vec(:,ind_mode,3) ]; % Cartesian vector representation      
    
    ax3 = MakeSphereAxis(p_res,t_res);     % Generate surface plot template
    ax3.CData = Grid2Mesh(vecnorm(Phi_Cart,2,2),p_res,t_res); % Plot vector norm
        
    % Add vector field directions
    hold on
    quiver3(e_r(:,1), e_r(:,2), e_r(:,3), Phi_Cart(:,1), Phi_Cart(:,2), Phi_Cart(:,3),...
            'Color','r','Linewidth',1.5,'AlignVertexCenters','on','AutoScaleFactor',1.2);
        
    % Add plot title
    title(['$\Phi_{i}^{(',num2str(l),',',num2str(m),')}$'],...
           'Interpreter','latex','Fontsize',30);    
   
        
    %%%%--- EXAMPLE 5: Tensor harmonic: (Psi_lm)_ij for chosen mode ---%%%%
    Psi_tt_lm_sel = PSI_tt_lm(:,ind_mode); % (theta,theta) component
    Psi_pt_lm_sel = PSI_pt_lm(:,ind_mode); % (phi,theta) component
    
    % Covariant determinant
    detPsi_lm = -(Psi_tt_lm_sel.^2 + PHI_tt_lm(:,ind_mode).^2); 
    
    % In-plane orientation of nematic directors
    [nemdir_psi_sc,nemdir_psi] = ...
        NemTensDir(Psi_tt_lm_sel,Psi_pt_lm_sel,detPsi_lm,e_r,e_t,e_p,nem_sc); 
    
    ax4 = MakeSphereAxis(p_res,t_res); % Generate surface plot template
    ax4.CData = Grid2Mesh(sqrt(abs(detPsi_lm)),p_res,t_res); % Colorcode with tensor determinant
    
    % Add nematic directors to plot
    hold on
    quiver3(nemdir_psi_sc(:,1), nemdir_psi_sc(:,2), nemdir_psi_sc(:,3),...
            nemdir_psi(:,1), nemdir_psi(:,2), nemdir_psi(:,3),...
            'Color','r','Linewidth',2,'ShowArrowHead','off');
    
    % Add plot title
    title(['$\Psi_{ij}^{(',num2str(l),',',num2str(m),')}$'],...
           'Interpreter','latex','Fontsize',30);
       
    
    %%%%--- EXAMPLE 6: Tensor harmonic: (Phi_lm)_ij for chosen mode ---%%%%
    Phi_tt_lm_sel = PHI_tt_lm(:,ind_mode); % (theta,theta) component
    Phi_pt_lm_sel = PHI_pt_lm(:,ind_mode); % (phi,theta) component
    
    % Covariant determinant (same for both tensors)
    detPhi_lm = detPsi_lm; 
    
    % In-plane orientation of nematic directors
    [nemdir_phi_sc,nemdir_phi] = ...
        NemTensDir(Phi_tt_lm_sel,Phi_pt_lm_sel,detPhi_lm,e_r,e_t,e_p,nem_sc); 
    
    ax5 = MakeSphereAxis(p_res,t_res); % Generate surface plot template
    ax5.CData = Grid2Mesh(sqrt(abs(detPhi_lm)),p_res,t_res); % Colorcode with tensor determinant
    
    % Add nematic directors to plot
    hold on
    quiver3(nemdir_phi_sc(:,1), nemdir_phi_sc(:,2), nemdir_phi_sc(:,3),...
            nemdir_phi(:,1), nemdir_phi(:,2), nemdir_phi(:,3),...
            'Color','r','Linewidth',2,'ShowArrowHead','off');
        
    % Add plot title
    title(['$\Phi_{ij}^{(',num2str(l),',',num2str(m),')}$'],...
           'Interpreter','latex','Fontsize',30);
    
end