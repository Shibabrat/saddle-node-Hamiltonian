function [kinetic_energy,isterminal,direction] = turning_point_event(t,x) 
    
    if abs(t) > 1
        isterminal = 1; % Halt integration 
    else
        isterminal = 0; % don't terminate within a short time
    end
    
    % The value that we want to be zero
    kinetic_energy = 0.5*(x(3).^2 + x(4).^2); 
    
    % The zero can be approached from either direction
    direction = 0; %0: all directions of crossing

     

end