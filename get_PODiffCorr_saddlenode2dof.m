function [x0po,t1] = get_PODiffCorr_saddlenode2dof(x0, par)

% [x0po,t1] = get_PODiffCorr_saddlenode2dof(x0, par)
% 
% Differential correction for the guess of the unstable periodic orbit initial
% condition. It keeps the initial x-position constant and varies the
% y-position value. The x-velocity is used for terminating the orbit at the
% half-period (may not work if the UPO is not symmetric), and y-velocity is
% used for testing convergence.
% 
% Shibabrat Naik (modified on 02-May-2019)
% parameters = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];
% 
% 


show = 1; % show = 1 to plot successive approximations (default=0)
% axesFontName = 'factory';
% axesFontName = 'Times New Roman';
% label_fs = 20; axis_fs = 30; % fontsize for publications 
label_fs = 10; axis_fs = 15; % fontsize for publications 
% set(0,'Defaulttextinterpreter','latex', ...
%     'DefaultAxesFontName', axesFontName, ...
%     'DefaultTextFontName', axesFontName, ...
%     'DefaultAxesFontSize', axis_fs, ...
%     'DefaultTextFontSize', label_fs, ...
%     'Defaultuicontrolfontweight','normal', ...
%     'Defaulttextfontweight','normal', ...
%     'Defaultaxesfontweight','normal');


% tolerances for integration and perpendicular crossing of x-axis
% MAXdxdot1 = 1.e-8 ; RelTol = 3.e-10; AbsTol = 1.e-10; 
MAXdydot1 = 1.e-12 ; RelTol = 3.e-14; AbsTol = 1.e-14; 
MAXdxdot1 = 1.e-12 ; RelTol = 3.e-14; AbsTol = 1.e-14; 

MAXattempt = 10;     	% maximum number of attempts before error is declared

dxdot1 	   = 1;         % to start while loop
dydot1 	   = 1;         % to start while loop
attempt    = 0;         % begin counting number of attempts
% y0(attempt) = 1;

while abs(dydot1) > MAXdydot1
% while abs(dxdot1) > MAXdxdot1
    
	if attempt > MAXattempt
		ERROR = 'Maximum iterations exceeded' ;
		disp(ERROR) ;
		break
    end
    
    y0 = x0;
    % Find first half-period crossing event
    TSPAN = [0 20]; % allow sufficient time for the half-period crossing event         
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol, ...
        'Events',@half_period_event); 
%     [tt,xx,t1,xx1,i1] = ode113(MODEL,TSPAN,x0,OPTIONS);
    [tt,xx,t1,xx1,i1] = ode113(@(t, y)saddlenode2dof(t, y, par), ...
                                TSPAN,x0,OPTIONS);
    
    
	x1 = xx1(end,1); 
    y1 = xx1(end,2); 
	dxdot1 = xx1(end,3); 
    dydot1  = xx1(end,4); 
%     plot3(xx(:,1),xx(:,2),xx(:,4),'-r');hold on
%     plot3(x1,y1,dydot1,'xb');
   
    
    % Compute the state transition matrix from the initial state to
	% the final state at the half-period event crossing
      
    % Events option not necessary anymore
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 
	[x,t,phi_t1,PHI] = stateTransitMat_saddlenode2dof(x0,t1(end),OPTIONS,par) ;

	attempt = attempt+1 ;
    
% 	ATTEMPT = sprintf('::poDifCor : iteration %d',attempt) ;
% 	disp(ATTEMPT) ;
     
    if show == 1
        e = get_TE_saddlenode2dof(x0,par);
        
%         plot3(x(:,1),x(:,2),x(:,4),'.-');
        plot3(x(:,1),x(:,2),x(:,4),'.-',x(:,1),x(:,2),-x(:,4),'.-'); 
        hold on;
        m = length(x) ;
        plot3(x(1,1),x(1,2),x(1,4),'r*');
        plot3(x(m,1),x(m,2),x(m,4),'go');
        
        
%         set(gca,'fontsize',18)
%         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
%         axis equal % destroyes the scale for successive plots
%         set(gca,'fontsize',label_fs)
%         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
%         axis equal % destroyes the scale for successive plots
%         xlabel('$q$','interpreter','latex','fontsize',axis_fs);
%         ylabel('$x$','interpreter','latex','fontsize',axis_fs);
%         zlabel('$p_x$','interpreter','latex','fontsize',axis_fs);
        xlabel('$x$','interpreter','latex','fontsize',axis_fs);
        ylabel('$y$','interpreter','latex','fontsize',axis_fs);
        zlabel('$p_y$','interpreter','latex','fontsize',axis_fs);
        
        title(['$\Delta E$ = ',num2str(mean(e) - 0)], ...
            'interpreter','latex','fontsize',axis_fs);
%         xlim([-15 15])
%         ylim([-15 15])
        pause(0.01) ;
        grid on
        box on

    end


%=========================================================================
% differential correction and convergence test, adjust according to
% the particular problem

    %compute acceleration values for correction term
    dVdx = ( (-2*sqrt(par(3))*x1 + par(4)*x1^2) - par(6)*(y1 - x1) );
    dVdy = ( par(5)^2*y1 + par(6)*y1 - par(6)*x1 );
    
    vxdot1 = -dVdx;
    vydot1 = -dVdy;
    
    %correction to the initial x0
%     y0(attempt) = dxdot1/(phi_t1(3,1) - phi_t1(4,1)*(vxdot1/vydot1));	
%     x0(1) = x0(1) - y0(attempt);
    
    %correction to the initial y0
    y0(attempt) = dydot1/(phi_t1(4,2) - phi_t1(3,2)*vydot1*(1/vxdot1)); 
	x0(2) = x0(2) - y0(attempt);
  
   
end

x0po = x0;
t1 = t1(end);

end