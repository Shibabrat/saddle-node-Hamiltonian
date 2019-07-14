function [eqPt] = get_eq_pts_saddlenode2dof(eqNum, parameters)

% get_eq_pts_saddlenode2dof solves for the equilibrium point assuming the
% system has KE + PE form and KE only depends on momenta.
% parameters = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];
    
    %fix the equilibrium point numbering convention here and make a
    %starting guess at the solution
    if 	eqNum == 1 
        x0 = [0; 0; 0; 0];  % EQNUM = 1, saddle  
    elseif 	eqNum == 2 
        x0 = [0.05; 0.2; 0; 0];  % EQNUM = 2, stable
    end
    
    %F(xEq) = 0 at the equilibrium point, solve using in-built function
    
%    options = optimoptions('fsolve','Display','iter'); % Option to display output
%    [eqPt,fval] = fsolve(@func_vec_field_eq_pt,x0,options) % Call solver
    
    [eqPt, fval, ~] = ....
        fsolve(@(x)func_vec_field_eq_pt(x,parameters), x0, ...
            optimset("jacobian","off")); % Call solver
 
end
function F = func_vec_field_eq_pt(x,par)
    

    dVdx = ( (-2*sqrt(par(3))*x(1) + par(4)*x(1)^2) - par(6)*(x(2) - x(1)) );
    dVdy = ( par(5)^2*x(2) + par(6)*x(2) - par(6)*x(1) );
    
    F = [   x(3)/par(1);
            x(4)/par(2);
            -dVdx;
            -dVdy];
    
end



