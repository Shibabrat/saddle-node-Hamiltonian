% Making a movie of the variation of the PES with the varying coupling
% parameter, EPSILON. 
% 
% Dependencies:
%   draw_SaddleNode_limitCurve(slice,fix_val,H0,mu,alpha,omega,epsi,col_curve)


MU = 0.1;
ALPHA = 2.0;
OMEGA = 1.0;
% EPSILON = 0.1;
EPSILON_c = (2*sqrt(MU)*OMEGA^2)/(OMEGA^2 - 2*sqrt(MU));

set(0,'DefaultLegendAutoUpdate','off')

frames = 60;
epsilon_vals = linspace(0,3,frames+1);
for k = 1:frames+1
    
%     draw_SaddleNode_limitCurve('qdot',0,-0.05,0.25,1,1.25,epsilon_vals(k),'g');

    EPSILON = epsilon_vals(k);
    % Evaluate the equilibrium point
    eqNum_2(1) = (2*sqrt(MU))/ALPHA - ...
        (OMEGA^2*epsilon_vals(k))/(ALPHA*(OMEGA^2 + EPSILON));
    eqNum_2(2) = (EPSILON/(OMEGA^2 + EPSILON))*eqNum_2(1);
    
    p1 = draw_SaddleNode_limitCurve('qdot',0,-0.05,0.25,1,1.25,epsilon_vals(k),'b');
    hold on
    p2 = draw_SaddleNode_limitCurve('qdot',0,-0.2,0.25,1,1.25,epsilon_vals(k),'g');
    p3 = draw_SaddleNode_limitCurve('qdot',0,0.0,0.25,1,1.25,epsilon_vals(k),'k');
    p4 = draw_SaddleNode_limitCurve('qdot',0,0.1,0.25,1,1.25,epsilon_vals(k),'r');
%     leg_limitcurve = legend([p1 p2 p3 p4], '$\mathcal{H} = -0.05$', ...
%         '$\mathcal{H} = -0.2$', '$\mathcal{H} = 0.0$', ...
%         '$\mathcal{H} = 0.1$');
    leg_limitcurve = legend([p1 p2 p3 p4], '$H_0 = -0.05$', ...
        '$H_0 = -0.2$', '$H_0 = 0.0$', '$H_0 = 0.1$');

    
%     leg1 = legend('$\bar{x}$','$\tilde{x}$','$\hat{x}$');
    set(leg_limitcurve,'Interpreter','latex');
    set(leg_limitcurve,'FontSize',15);
    
%     if epsilon_vals(k) < EPSILON_c
%         scatter(0,0,40,'+r');
%     end
%     scatter(eqNum_2(1),eqNum_2(2),20,'or','filled');
    
    title(['$\epsilon = $',num2str(epsilon_vals(k))])
    
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.02];
    axis equal

    fig.PaperPositionMode = 'auto';
%     print(['./snapshots_movie/output_',num2str(k),'.png'], ...
%         '-dpng','-r600')
    print(['../../data/saddlenode-ham/2dof_movie_v2/epsilon_', ...
        num2str(k),'.png'], '-dpng','-r600')
   
    close all
end






%%
















