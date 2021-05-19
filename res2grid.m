function grid = res2grid(pres,tres)
% Creates linear list of spherical coordinate grid in the format
% needed to evaluate harmonics at these given points
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alexander Mietke, 05/19/2021
%   alexanmie@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [T, P] = meshgrid(tres:tres:180,0:pres:360);
    
    % Convert to radian
    T = 2 * pi * T / 360;
    P = 2 * pi * P / 360;
    
    if mod(360,pres)==0
        % Remove doubled phi coordinates
        P(end,:) = [];
        T(end,:) = [];
    end
    
    if mod(180,tres)==0
        % Remove repeated pole coordinates
        T(:,end) = [];
        P(:,end) = [];
    end
    
    % Generate grid coordinates including poles
    grid = [0, 0; P(:), T(:); 0, pi];
    
    
end