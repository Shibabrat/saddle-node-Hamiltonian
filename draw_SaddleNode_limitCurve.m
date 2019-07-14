function plt_hand = draw_SaddleNode_limitCurve(slice,fix_val,H0,mu,alpha,omega,epsi,col_curve)   
    
% @author: Víctor J. García-Garrido
% 
% Modified by Shibabrat Naik on 30-Apr-2019

    axes_fontsize = 18;

    % remember that qdot = p / mass ; xdot = px / mass  
    if strcmp(slice,'qdot')
        % qdot = fix_val ; xdot = 0
        fix_val = 0;       
        f = @(q,x) -sqrt(mu).*q.^2 + (alpha./3).*q.^3 + ...
            (1/2).*omega.^2.*x.^2 + (epsi./2).*(x-q).^2 - H0;
        
        q0 = -1.5;
        q1 = 2.5;
        x0 = -1.25;
        x1 = 1.25;
        
        plt_hand = fimplicit(f,[q0 q1 x0 x1],'Color',col_curve,'LineWidth',2);
        
        axis([q0 q1 x0 x1]);
        xlabel('$q$','FontSize',axes_fontsize,'Interpreter','Latex');
        ylabel('$x$','FontSize',axes_fontsize,'Interpreter','Latex','Rotation',90);
        % dist_x = x1 - x0;
        % dist_y = y1 - y0;
        % daspect([1 dist_y/dist_x 1]);      
        
    elseif strcmp(slice,'q')
        % q = fix_val ; qdot = 0
        q_well = (2*sqrt(mu))/alpha - ((omega^2 * epsi)/(alpha*(omega^2+epsi)));
        q0 = q_well;
        f = @(x,px) (1./2).*px.^2 - sqrt(mu).*q0.^2 + (alpha./3).*q0.^3 + ...
            (1/2).*omega.^2.*x.^2 + (epsi./2).*(x-q0).^2 - H0;
      
        [x0,x1,px0,px1] = x_px_vals_qwell(H0,mu,alpha,omega,epsi);
        
        x0 = x0 - 0.02;
        x1 = x1 + 0.02;
        px0 = px0 - 0.02;
        px1 = px1 + 0.02;
        
        fimplicit(f,[x0 x1 px0 px1],'Color',col_curve,'LineWidth',2)      
        xlabel('$x$','FontSize', axes_fontsize,'Interpreter','Latex');
        ylabel('$p_x$','FontSize', axes_fontsize,'Interpreter','Latex','Rotation',0);
        axis([x0 x1 px0 px1]);

    elseif strcmp(slice,'x')
        % x = fix_val ; xdot = 0
        x0 = fix_val;
        
        f = @(q,p) (1./2).*p.^2 - sqrt(mu).*q.^2 + (alpha./3).*q.^3 + ...
            (1/2).*omega.^2.*x0.^2 + (epsi./2).*(x0-q).^2 - H0;
        
        q0 = -0.05;
        q1 = 2.25;
        p0 = -1.6;
        p1 = 1.6;
            
        fimplicit(f,[q0 q1 p0 p1],'Color',col_curve,'LineWidth',2)
        xlabel('$q$','FontSize',axes_fontsize,'Interpreter','Latex');
        ylabel('$p$','FontSize',axes_fontsize,'Interpreter','Latex','Rotation',0);
        axis([q0 q1 p0 p1]);
%         dist_y = y1 - y0;
%         dist_py = py1 - py0;
%         daspect([1 dist_py/dist_y 1]);  
    else
        % x = 0 ; q = 0
        f = @(p,px) (1 ./ 2) .* (px.^2 + p.^2) - H0;
        epsi = 0.05;
        px0 = - sqrt(2*H0) - epsi;
        px1 = sqrt(2*H0) + epsi;
        p0 = - sqrt(2*H0) - epsi;
        p1 = sqrt(2*H0) + epsi;
        fimplicit(f,[p0 p1 px0 px1],'Color',col_curve,'LineWidth',2)
        xlabel('$p$','FontSize',axes_fontsize,'Interpreter','Latex');
        ylabel('$p_x$','FontSize',axes_fontsize,'Interpreter','Latex','Rotation',0);
        axis equal
    end
    
%     set(gca,'FontSize',axes_fontsize);
     
    xl = get(gca,'XLabel');
    xlFontSize = get(xl,'FontSize');
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', axes_fontsize);
    set(xl, 'FontSize', xlFontSize + 15);
    
     
    yl = get(gca,'YLabel');
    ylFontSize = get(yl,'FontSize');
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', axes_fontsize);
    set(yl, 'FontSize', ylFontSize + 15);
    
    
%     xlim([-1.5 1.5]);
%     ylim([-1.5 1.5])
    
end





