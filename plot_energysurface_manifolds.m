%%

label_fs = 30; % fontsize for publications
eqNum = 1;
MASS_A = 1.0; MASS_B = 1.0; % Mass-weighted momenta
MU = 4.0;
OMEGA = 3.0;
ALPHA = 1.0;
EPSILON = 5.0;
deltaE = 0.5;
data_path = ...
    'data-figures/manifolds/coupled_mu4_omega3_epsilon5_deltaEnergy5e-1/manifold-energysurface/';
% data_path = ...
%     'data-figures/manifolds/uncoupled_mu4_omega3_epsilon5_deltaEnergy5e-1/manifold-energysurface/';

parameters = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];


upo_lw = 3;

po_target_file = ['upo_','x0po_T_energyPO_eqPt',num2str(eqNum), ...
                    '_DelE',num2str(deltaE),'_saddlenode2dof.txt'];
                
xx = importdata([data_path, po_target_file]);

plot3([xx(:,1);flip(xx(:,1))], [xx(:,2);flip(xx(:,2))], ...
        [xx(:,4);flip(-xx(:,4))],'-k','Linewidth',upo_lw);
grid on
box on
hold on


fs = plot_energysurface(deltaE,parameters,0.3);

% fimplicit3(@(x,y,z) -sqrt(parameters(3))*x.^2 + (parameters(4)/3)*x.^3 ...
%             + 0.5*parameters(5)^2*y.^2 + 0.5*parameters(6)*(y - x).^2 ...
%             - (0 + deltaE), 'm', [-2 2 -2 2 -4 -3.95], ...
%             'MeshDensity', 100)

% view(62,43)
view(15,20)
xlabel('$x$','interpreter','latex','fontsize',label_fs);
ylabel('$y$','interpreter','latex','fontsize',label_fs);
zlabel('$p_y$','interpreter','latex','fontsize',label_fs);
% print('-dpng', '-r600','temp.png')


%%