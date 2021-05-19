function [NemDir_sc, NemDir] = NemTensDir(Qtt,Qpt,detQ,e_r,e_t,e_p,nem_sc)
% Calculates local in-plane angle of the nematic director field
% encoded by the trace-free symmetric harmonic tensor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine local in-plane angle
    psi_dummy = Qtt./sqrt(abs(detQ)); 
    nem_angle_psi = 0.5*sign(Qpt).*acos(psi_dummy);
    NemAngle = nem_angle_psi + 0.5*pi*(Qtt<0).*(Qpt==0);
    
    % Pole correction for Psi tensors
    NemAngle = NemAngle + ...
        0.25*pi*(Qtt==0).*(Qpt==0).*[1; zeros([length(Qtt)-2,1]);-1]; 
    
    % Compute local tangent vector accordingly
    NemDir = cos(repmat(NemAngle,[1,3])).*e_t + sin(repmat(NemAngle,[1,3])).*e_p; 
    
    % Center nematic directors for better visualization
    nem_start = e_r - nem_sc*NemDir;
    nem_end = e_r + nem_sc*NemDir;
    nem_start(isnan(nem_start)) = [];
    nem_end(isnan(nem_end)) = [];
    nr_nonan = length(nem_start(:))/3;
    NemDir_st = reshape(nem_start,[nr_nonan,3]);
    NemDir_end = reshape(nem_end,[nr_nonan,3]);
    
    % Prepare information to be used by quiver plot
    NemDir_sc = NemDir_st; % Quiver symbols plotted at start coordinate
    NemDir = NemDir_end - NemDir_sc; % Nematic director
end