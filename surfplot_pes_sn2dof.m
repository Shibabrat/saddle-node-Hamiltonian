function surfplot_pes_sn2dof(numpts)

    set(0, ...
    'Defaulttextinterpreter','latex', ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesFontSize', 20, ...
    'DefaultTextFontSize', 20, ...
    'DefaultAxesLabelFontSize', 2.0, ...
    'Defaultuicontrolfontweight','default', ...
    'Defaulttextfontweight','default', ...
    'Defaultaxesfontweight','default');


    xi = -10; xf = 10;
    yi = -10; yf = 10;
    xGrid = linspace(xi, xf, numpts);
    yGrid = linspace(yi, yf, numpts);
    [xMesh, yMesh] = meshgrid(xGrid, yGrid);

    % parameters = [MASS_A MASS_B MU ALPHA OMEGA EPSILON]
%     parameters = [1 1 4 1 3 0];
    parameters = [1 1 4 1 3 5];
    peMesh = get_peMesh(xMesh, yMesh, parameters);
    
    filename2save = ['pes2dof_mu',num2str(parameters(3)), ...
        '_alpha',num2str(parameters(4)), ...
        '_omega',num2str(parameters(5)), ...
        '_epsilon',num2str(parameters(6))];
    
    figure1 = figure();
    surfc(xMesh, yMesh, peMesh)
%     colormap winter
    shading interp
    view(-25,25)
    xlabel('$x$')
    ylabel('$y$')
%     legend({'$\mathcal{D} = 10.67 \, , \; \mathcal{F} = 30.09$'}, ...
%         'Interpreter','latex','Location','NorthWest')
    legend({'$\mathcal{D} = 0.08 \, , \; \mathcal{F} = 59.71$'}, ...
        'Interpreter','latex','Location','NorthWest')
    print(figure1, '-dpng','-r300',filename2save)

    
%     parameters = [1 1 4 2 3 0];
    parameters = [1 1 4 2 3 5];
    peMesh = get_peMesh(xMesh, yMesh, parameters);
    
    filename2save = ['pes2dof_mu',num2str(parameters(3)), ...
        '_alpha',num2str(parameters(4)), ...
        '_omega',num2str(parameters(5)), ...
        '_epsilon',num2str(parameters(6))];
    
    figure2 = figure();
    surfc(xMesh, yMesh, peMesh)
%     colormap winter
    shading interp
    view(-25,25)
    xlabel('$x$')
    ylabel('$y$')
%     legend({'$\mathcal{D} = 2.67 \, , \; \mathcal{F} = 55.20$'}, ...
%         'Interpreter','latex','Location','NorthWest')
    legend({'$\mathcal{D} = 0.02 \, , \; \mathcal{F} = 85.56$'}, ...
        'Interpreter','latex','Location','NorthWest')
    print(figure2, '-dpng','-r300',filename2save)
    
%     parameters = [1 1 4 5 3 0];
    parameters = [1 1 4 5 3 5];
    peMesh = get_peMesh(xMesh, yMesh, parameters);
    
    filename2save = ['pes2dof_mu',num2str(parameters(3)), ...
        '_alpha',num2str(parameters(4)), ...
        '_omega',num2str(parameters(5)), ...
        '_epsilon',num2str(parameters(6))];
    
    figure3 = figure();
    surfc(xMesh, yMesh, peMesh)
%     colormap winter
    shading interp
    view(-25,25)
    xlabel('$x$')
    ylabel('$y$')
%     legend({'$\mathcal{D} = 0.43 \, , \; \mathcal{F} = 142.06$'}, ...
%         'Interpreter','latex','Location','NorthWest')
    legend({'$\mathcal{D} = 0.00 \, , \; \mathcal{F} = 171.27$'}, ...
            'Interpreter','latex','Location','NorthWest')
    print(figure3, '-dpng','-r300',filename2save)
    
end

function peMesh = get_peMesh(xMesh, yMesh, parameters)

    peMesh = -sqrt(parameters(3))*xMesh.^2 + (parameters(4)/3)*xMesh.^3 + ...
            0.5*parameters(5)^2*yMesh.^2 + 0.5*parameters(6)*(yMesh - xMesh).^2;
    
    
end
