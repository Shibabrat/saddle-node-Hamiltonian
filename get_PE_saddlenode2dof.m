function pot_energy = get_PE_saddlenode2dof(q,par)

%    par = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];
    
    pot_energy = -sqrt(par(3))*q(1)^2 + (par(4)/3)*q(1)^3 + ...
        0.5*par(5)^2*q(2)^2 + 0.5*par(6)*(q(2) - q(1))^2;
                
    
end