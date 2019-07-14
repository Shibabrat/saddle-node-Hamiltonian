function [q,x] = SNode_wellLocation(mu,alpha,omega,epsi)
    
    q = (2.*sqrt(mu)./alpha) - (epsi.*omega.^2)./(alpha.*(epsi + omega.^2));
    x = (epsi./(epsi + omega.^2)) .* q;

end

