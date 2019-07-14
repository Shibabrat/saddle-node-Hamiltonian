function Df = get_Jacobian(xEq)

    global OMEGA_X OMEGA_Y DELTA
    
    syms x y vx vy 
    
    xe = xEq(1);
    ye = xEq(2);
    vxe = xEq(3);
    vye = xEq(4);
    
    %Use vector differentiation 
    f(x,y,vx,vy) = [vx; vy; -(OMEGA_X^2*x + DELTA*y^2); 
        -(OMEGA_Y^2*y + 2*DELTA*x*y)];
    DfMat(x,y,vx,vy) = jacobian(f,[x y vx vy]);
    Df = double(DfMat(xe,ye,vxe,vye));
    
end




