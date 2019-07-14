function [xDot,isterminal,direction] = full_period_event(t,x) 

    global y0
    
    % The value that we want to be zero
%     xDot = x(3);
%     size(x)
    % dDSQdt is the derivative of the equation for current distance. Local
    % minimum/maximum occurs when this value is zero.
    dDSQdt = 2 * ((x(1:2)-y0(1:2))' * x(3:4)); 
    xDot = [dDSQdt; dDSQdt]; 
    
    if abs(t) > 1e-2 
        isterminal = [0; 1]; % Halt integration at the maximum
    else
        isterminal = [0; 0]; % don't terminate within a short time
    end
    
    % The zero can be approached from either direction
%     direction = 0; %0: all directions of crossing
    direction  = [1; -1];          % [local minimum, local maximum]
    
    
     

end