% Parameters for the 1 DoF model
MU = 6;
ALPHA = 4;
filename2save = ['output_MU_',num2str(MU),'_ALPHA_',num2str(ALPHA)];

% Parameters of the fimplicit plot
% limits = [-4*sqrt(MU)/ALPHA 4*sqrt(MU)/ALPHA];
limits = [-8 8];
mesh_density = 201;
line_width = 3;
energy_reaction = 2;

set(0, ...
'Defaulttextinterpreter','latex', ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesFontSize', 25, ...
'DefaultTextFontSize', 25, ...
'DefaultAxesLabelFontSize', 2.0, ...
'Defaultuicontrolfontweight','default', ...
'Defaulttextfontweight','default', ...
'Defaultaxesfontweight','default');

close all
figure1 = figure('Position', [50 50 700 700]);

fimplicit(@(q,p)(0.5*p.^2 - MU^(0.5).*q.^2 + ALPHA*(q.^3/3) - 0), limits, ...
            'MeshDensity',mesh_density, ...
            'LineWidth',line_width,'Color','black')
hold on
react_traj = fimplicit(@(q,p)(0.5*p.^2 - MU^(0.5).*q.^2 + ALPHA*(q.^3/3) - ...
                            energy_reaction), limits, ...
                            'MeshDensity',mesh_density, ...
                            'LineWidth',line_width,'Color','red');
        
fimplicit(@(q,p)(0.5*p.^2 - MU^(0.5).*q.^2 + ALPHA*(q.^3/3) + 5), limits, ...
            'MeshDensity',mesh_density, ...
            'LineWidth',line_width,'Color','blue')
        
nonreact_traj = fimplicit(@(q,p)(0.5*p.^2 - MU^(0.5).*q.^2 + ALPHA*(q.^3/3) + 1), limits, ...
                            'MeshDensity',mesh_density, ...
                            'LineWidth',line_width,'Color','green');

                        
% Dividing surface and NHIM
nhim = scatter(0,0,400,'m','Marker','+','Linewidth',5);

ds = scatter([0,0],[-sqrt(2*energy_reaction),sqrt(2*energy_reaction)],100,'c', ...
            'Marker','o', ...
            'MarkerFaceColor', 'g', ...
            'Linewidth',5);

% leg_hand = legend([react_traj nonreact_traj nhim ds],'Reactive','Non-reactive', ...
%             'NHIM', 'DS', 'Location', [0.6, 0.78, 0.3, 0.1]);
% set(leg_hand,'FontSize',20);

ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.015 0.02];
axis equal


% annotation(figure1,'textarrow',[0.615, 0.52],[0.78875, 0.53], ...
%             'Color','m','String','NHIM');

% xlabel('$q$')
% ylabel('$p$','Rotation',0)
% print(figure1, '-dpng','-r600','output.png')

fig.PaperPositionMode = 'auto';
% print('output','-dpng','-r0')
print('-bestfit',filename2save,'-dpdf')

%%




