function M = Grid2Mesh(f,pres,tres)
% Rearrange linear vector of function values on the sphere 
% into a meshgrid-type matrix needed for surface plotting
% The reshaping is seperately performed on each colomnu of f
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of phi- and theta-grid points
    n_phi = 360/pres + 1;
    n_theta = 180/tres + 1;
    
    % Empty output matrix
    M = zeros([n_theta,n_phi,size(f,2)]);
    
    for i = 1:size(f,2) % Loop through columns of f
        f_curr = f(:,i);
        
        % Rearrange non-pole points and add double-row at phi = 0, 2*pi
        M_np = reshape(f_curr(2:(end-1)),[n_phi-1,n_theta-2]);
        M_np = [M_np; M_np(1,:)];
        M_np = M_np';

        % Create meshgrid format array with function data
        M(:,:,i) = [ repmat(f_curr(1),[1,n_phi]); ...        % pole theta = 0
                     M_np;  ...                              % 0 < theta < pi
                     repmat(f_curr(end),[1,n_phi]) ];        % % pole theta = pi
    end
end