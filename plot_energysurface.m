function fs = draw_energysurf_saddlenode2dof(H_val,par,alpha)

    % plot properties
    axesFontName = 'factory';
    % axesFontName = 'Times New Roman';
    axFont = 18;
    textFont = 18;
    labelFont = 25;
    lw = 2;    
    set(0,'Defaulttextinterpreter','latex', ...
        'DefaultAxesFontName', axesFontName, ...
        'DefaultTextFontName', axesFontName, ...
        'DefaultAxesFontSize',axFont, ...
        'DefaultTextFontSize',textFont, ...
        'Defaultuicontrolfontweight','normal', ...
        'Defaulttextfontweight','normal', ...
        'Defaultaxesfontweight','normal');
    
    esurf = @(x,y,px) 2*((1/(2*par(1)))*px.^2 + ...
        -sqrt(par(3))*x.^2 + (par(4)/3)*x.^3 + ...
        0.5*par(5)^2*y.^2 + 0.5*par(6)*(y - x).^2 - H_val);
    
    
%     rgb_col = [0/255 255/255 0/255];
    lightGrey1   = [0.85 0.85 0.85];
    lightGrey2   = [0.7 0.7 0.7];

    darkGrey1  = [0.4 0.4 0.4];
    darkGrey2  = [0.2 0.2 0.2];
    rgb_col = darkGrey1;
     
    xi = -2; xf = 2;
    yi = -2; yf = 2;
    pxi = -2; pxf = 2;
%     pxi = -1; pxf = 0.5;
    
    xi = -2; xf = 10;
    yi = -2; yf = 10;
    pxi = -5; pxf = 5;


%     xi = -1; xf = 1;
%     yi = -1; yf = 1;
%     pxi = -1; pxf = 1;

    fs = fimplicit3(esurf,[xi xf yi yf pxi pxf],...
        'EdgeColor','none','MeshDensity',100,'FaceAlpha',alpha,...
        'FaceColor',rgb_col);
  
%     xlabel('$x$','FontSize',labelFont,'Interpreter','Latex');
%     ylabel('$y$','FontSize',labelFont,'Interpreter','Latex');
%     zlabel('$p_y$','FontSize',labelFont,'Interpreter','Latex');
    light;
    
    % fs.FaceAlpha = alpha;
    fs.AmbientStrength = 0.5;
    
end

