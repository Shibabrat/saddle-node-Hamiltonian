function e = get_TE_saddlenode2dof(orbit, parameters)

%   get_total_energy_saddlenode2dof computes the total energy of an input orbit
%   (represented as M x N with M time steps and N = 4, dimension of phase
%   space for the model) for the 2 DoF saddle node Hamiltonian.
% 
%   Orbit can be different initial conditions for the periodic orbit of
%   different energy. When trajectory is input, the output energy is mean.
%   parameters = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];

    
    x = orbit(:,1);
    y = orbit(:,2);
    px = orbit(:,3);
    py = orbit(:,4);
    
    e = (1/(2*parameters(1)))*px.^2 + (1/(2*parameters(2)))*py.^2 + ...
            get_PE_saddlenode2dof([x, y],parameters);
    
%     if length(e) > 1 
%         e = mean(e);
%     end
        
    
end
