function [xDot,isterminal,direction] = intersect_sosatxw_event(t,x) 

% [xDot,isterminal,direction] = intersect_sosatxw_event(t,x)
% 
% Function to catch event of crossing x = x_w crossing of a trajectory
% parameters = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];

    global parameters
    ALPHA = parameters(4);
    MU = parameters(3);
    OMEGA = parameters(5);
    EPSILON = parameters(6);
    
    x_w = (2*sqrt(MU))/ALPHA - ...
        (OMEGA^2*EPSILON)/(ALPHA*(OMEGA^2 + EPSILON));
    
    xDot = x(1) - x_w;
        
    if abs(t) > 1e-2 
        isterminal = 1; % Halt integration 
    else
        isterminal = 0; % don't terminate within a short time
    end
    
    % The zero can be approached from either direction
%     direction = 0; %0: all directions of crossing 
%     direction = -1; %works for x(3) < 0 intersection of unstable manifolds
%     direction = 1; %works for x(3) > 0 intersection of unstable manifolds
%     direction = 1; %works for x(3) < 0 intersection of stable manifolds
    direction = -1; %works for x(3) > 0 intersection of stable manifolds
    
end