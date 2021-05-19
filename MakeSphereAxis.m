function ax = MakeSphereAxis(p_res,t_res)
% Create a sphere surface plot template 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Create meshgrid for given angle resolution
    [P, T] = meshgrid(pi*(0:p_res:360)/180,pi*(0:t_res:180)/180);
    
    % Create new figure
    figure;

    % Create plot and set layout
    ax = surf(sin(T).*cos(P),sin(T).*sin(P),cos(T));       
    ax.FaceColor = 'interp';
    ax.LineStyle = 'none';
    ax.FaceLighting = 'none';
    colormap(parula(512));
    camlight;
    shading interp;
    axis equal;
    axis off;
end