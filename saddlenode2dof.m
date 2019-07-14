function xDot = saddlenode2dof(t, x, par)

    xDot = zeros(length(x),1);
    
    dVdx = ( (-2*sqrt(par(3))*x(1) + par(4)*x(1)^2) - par(6)*(x(2) - x(1)) );
    dVdy = ( par(5)^2*x(2) + par(6)*x(2) - par(6)*x(1) );
    
    xDot(1) = x(3)/par(1);
    xDot(2) = x(4)/par(2);
    xDot(3) = -dVdx;
    xDot(4) = -dVdy;

end