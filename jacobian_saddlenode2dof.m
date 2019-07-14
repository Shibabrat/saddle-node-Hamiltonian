function Df = jacobian_saddlenode2dof(eqPt, par)

     % par = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];

%     global OMEGA_X OMEGA_Y DELTA
    
    syms x y vx vy 
    
    xe = eqPt(1);
    ye = eqPt(2);
    vxe = eqPt(3);
    vye = eqPt(4);
    
    % Use vector differentiation 
    f(x,y,vx,vy) = [vx/par(1); 
                    vy/par(2); 
                   -( (-2*sqrt(par(3))*x + par(4)*x^2) - par(6)*(y - x) ); 
                   -( par(5)^2*y + par(6)*y - par(6)*x )];
    
    DfMat(x,y,vx,vy) = jacobian(f,[x y vx vy]);
    
    Df = double(DfMat(xe,ye,vxe,vye));
    
end





