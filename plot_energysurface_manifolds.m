%%
data_path = 'data-figures/manifolds/'

upo_linewidth = 4;


xx = importdata(['upo_',po_target_file]);


plot3([xx(:,1);flip(xx(:,1))], [xx(:,2);flip(xx(:,2))], ...
        [xx(:,4);flip(-xx(:,4))],'-k','Linewidth',upo_lw);
grid on
box on
hold on


fs = draw_energysurf_saddlenode2dof(deltaE,parameters,0.3);
view(62,43)



%%