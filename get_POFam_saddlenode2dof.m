function [x0po,T] = ...
    get_POFam_saddlenode2dof(eqNum,Ax1,Ax2,nFam,po_fam_file,parameters) 

% [x0po,T] =
% get_POFam_saddlenode2dof(eqNum,Ax1,Ax2,nFam,po_fam_file,parameters)
% 
% Generates a family of unstable periodic orbits (UPOs) given a pair of
% seed initial conditions from a data file, while targeting a specific
% periodic orbit.
% 
% Shibabrat Naik (modified on 02-May-2019)

    % delt = guessed change in period between successive orbits in family
    delt = -1.e-10 ;   % <==== may need to be changed, delt = -1.e-12 ;   

    N = 4 ; % dimension of phase space

    x0po = zeros(nFam,N) ;
    T    = zeros(nFam,1) ;
    energyPO = zeros(nFam,1) ;
    
    [x0poGuess1,TGuess1] = ...
            get_POGuessLinear_saddlenode2dof(eqNum,Ax1,parameters)
    [x0poGuess2,TGuess2] = ...
            get_POGuessLinear_saddlenode2dof(eqNum,Ax2,parameters)

    % Get the first two periodic orbit initial conditions
    iFam = 1 ;
    FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
    disp(FAMNUM) ;
    [x0po1,tfpo1] = get_PODiffCorr_POFam(x0poGuess1, parameters) ; 
    energyPO(iFam) = get_TE_saddlenode2dof(x0po1', parameters);


    iFam = 2 ;
    FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
    disp(FAMNUM) ;
    [x0po2,tfpo2] = get_PODiffCorr_POFam(x0poGuess2, parameters) ;               
    energyPO(iFam) = get_TE_saddlenode2dof(x0po2', parameters);


    x0po(1:2,1:N) = [x0po1(:)'  ; x0po2(:)' ];
    T   (1:2)     = [2*tfpo1    ; 2*tfpo2     ]; 

    %Generate the other members of the family using numerical continuation
    for iFam = 3:nFam

        FAMNUM = sprintf('::poFamGet : number %d', iFam) ;
        disp(FAMNUM) ;
        
        dx  = x0po(iFam-1,1) - x0po(iFam-2,1) ;
        dy  = x0po(iFam-1,2) - x0po(iFam-2,2) ;
%         dOrbit = x0po(iFam-1,:) - x0po(iFam-2,:);
        dt  = T(iFam-1) - T(iFam-2) ;

%         x0po_g = x0po(iFam-1,:)' + dOrbit';
        t1po_g =   (T(iFam-1) + dt)/2 + delt ;
        x0po_g = [ (x0po(iFam-1,1) + dx) (x0po(iFam-1,2) + dy) 0 0] ;
%         t1po_g = T(iFam-2) + abs(dt);

      % differential correction takes place in the following function
        [x0po_iFam,tfpo_iFam] = get_PODiffCorr_POFam(x0po_g, parameters) ; 
%         [x0po_iFam,tfpo_iFam] = ...
%             po_auto_shooting_ball_rolling(x0po_g', t1po_g);

        x0po(iFam,1:N) = x0po_iFam ;
%         T   (iFam)     = 2*t1_iFam ;	
        T(iFam)        = 2*tfpo_iFam;

        energyPO(iFam) = get_TE_saddlenode2dof(x0po(iFam,:), parameters) ;
        
        if mod(iFam,10) == 0
            dum = [x0po T] ;
    %     save x0po_T.dat -ascii -double dum
        end

    end

    dum = [x0po T energyPO] ;
    save(po_fam_file,'dum','-ascii','-double');

end
function [x0po,t1] = get_PODiffCorr_POFam(x0, par)

% global eSaddle
% global OMEGA_X OMEGA_Y DELTA

% set show = 1 to plot successive approximations (default=0)
show = 1 ;
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
MAXdxdot1 = 1.e-12 ; RelTol = 3.e-14; AbsTol = 1.e-14; 
MAXdydot1 = 1.e-12 ; RelTol = 3.e-14; AbsTol = 1.e-14; 

MAXattempt = 20;     	% maximum number of attempts before error is declared

dxdot1 	   = 1;         % to start while loop
dydot1 	   = 1;         % to start while loop
attempt    = 0;         % begin counting number of attempts
% y0(attempt) = 1;

% while abs(dxdot1) > MAXdxdot1
while abs(dydot1) > MAXdydot1

	if attempt > MAXattempt
		ERROR = 'Maximum iterations exceeded' ;
		disp(ERROR) ;
        break
    end
    
    y0 = x0;
    % Find first half-period crossing event
%     MODEL = 'deleonberne2dof';  
    TSPAN = [0 60];        % allow sufficient time for the half-period crossing event         
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol, ...
        'Events',@half_period_event); 
%     OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol, ...
%         'Events',@turning_point_event); 
    
%     [tt,xx,t1,xx1,i1] = ode113(MODEL,TSPAN,x0,OPTIONS);
    [tt,xx,t1,xx1,i1] = ode113(@(t, y)saddlenode2dof(t, y, par), ...
                                TSPAN,x0,OPTIONS);
                            
%     [tt,xx,t1,xx1,i1] = ode113(@(t, y)saddlenode2dof(t, y, par), ...
%                                 TSPAN,xx1(end,:),OPTIONS);
                     
	x1 = xx1(end,1); 
    y1 = xx1(end,2); 
	dxdot1 = xx1(end,3); 
    dydot1  = xx1(end,4); 
%     plot3(xx(:,1),xx(:,2),xx(:,3),'-r');hold on
%     plot3(x1,y1,dxdot1,'xb');
   
    
    % Compute the state transition matrix from the initial state to
	% the final state at the half-period event crossing
      
    % Events option not necessary anymore
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 
	[x,t,phi_t1,PHI] = stateTransitMat_saddlenode2dof(x0,t1(end),OPTIONS,par) ;

	attempt = attempt+1 ;
    
	ATTEMPT = sprintf('::poDifCor : iteration %d',attempt) ;
	disp(ATTEMPT) ;
     
    if show == 1
        e = get_TE_saddlenode2dof(x,par);
        plot3(x(:,1),x(:,2),x(:,4),'.-',x(:,1),x(:,2),-x(:,4),'.-'); 
        hold on;
        m = length(x) ;
        plot3(x(1,1),x(1,2),x(1,4),'g*');
        plot3(x(m,1),x(m,2),x(m,4),'ro');
        
        
%         set(gca,'fontsize',18)
%         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
%         axis equal % destroys the scale for successive plots
%         set(gca,'fontsize',label_fs)
%         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
%         axis equal % destroyes the scale for successive plots
        xlabel('$q$','interpreter','latex','fontsize',axis_fs);
        ylabel('$x$','interpreter','latex','fontsize',axis_fs);
        zlabel('$p_x$','interpreter','latex','fontsize',axis_fs);
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














