function [xDot,isterminal,direction] = half_period_event(t,x) 
    
    if abs(t) > 1e-2 
        isterminal = 1; % Halt integration 
    else
        isterminal = 0; % don't terminate within a short time
    end
    
    % The value that we want to be zero
    xDot = x(3); 
    
    % The zero can be approached from either direction
    direction = 0; %0: all directions of crossing

end