% Load the intersection of the manifolds with the SOS: q = q_w and p > 0 or
% p < 0
data_path = '../../data/saddlenode/2dof-set2/';
% data_path = './';

unmanipos_intersect_sos = importdata([data_path, ...
    'xeU1_unstable_branch1_eqPt1_DelE0.05_saddlenode2dof_trajs1000_p-pos.txt']);
smanipos_intersect_sos = importdata([data_path, ...
    'xeU1_stable_branch1_eqPt1_DelE0.05_saddlenode2dof_trajs1000_p-pos.txt']);

% Integrate the points on this SOS: backward for stable and forward for
% unstable manifold
% Setting up parameters and global variables
N = 4;          % dimension of phase space
MASS_A = 8.0; MASS_B = 8.0; % De Leon, Marston (1989)
EPSILON_S = 1.0;
D_X = 10.0;
% ALPHA = 2.30;
% LAMBDA = 1.95;
ALPHA = 1.00;
LAMBDA = 1.5;
par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];
total_energy = 1.01;

tb = -20;
tf = 0;
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14); % high accuracy
eventSwitch = 'on';
numTimesteps = 500;

iterate = {}; % cell array to store interate of the first intersection
curr_intersection = smanineg_intersect_sos;
max_iterate = 3;
j = 1;
while j <= max_iterate
    curr_iterate = [];
    for i = 1:size(curr_intersection,1)
        [x,t,te,xe,ie] = get_traj_sos_deleonberne(curr_intersection(i,:), ...
            tb, tf, OPTIONS, eventSwitch, numTimesteps, par);

        curr_iterate = [curr_iterate; xe(end,:)];
    end
    curr_intersection = curr_iterate;
    iterate{j} = curr_iterate;
    plot(iterate{j}(:,1),iterate{j}(:,3),'-')
    hold on
    j = j + 1;
end

%%

% domain = [-0.025 0.025 -0.5 0.5];
domain = [-0.25 0.25 -5 5];
fimplicit(@(x,px) par(4)*( 1 - exp(-par(5).*x) ).^2 - ...
            exp(-par(6)*par(5).*x) + par(3) + px.^2/(2*par(1)) - ...
            total_energy, domain, '-k','LineWidth',1,'MeshDensity', 1000)
    
        
% axis equal
daspect([(domain(2)-domain(1)) (domain(4)-domain(3)) 1])
xlabel('$x$')
ylabel('$p_x$')
xticks([domain(1) 0 domain(2)])
xticklabels({num2str(domain(1)),'0',num2str(domain(2))})
set(gca,'TickDir','out','TickLength',[0.02 0.02]); % The only other option is 'in'

hold on

% plot(unmanineg_intersect_sos(:,1), unmanineg_intersect_sos(:,3),'.m')
% plot(unmanipos_intersect_sos(:,1), unmanipos_intersect_sos(:,3),'.m')
plot(smanineg_intersect_sos(:,1), smanineg_intersect_sos(:,3),'-m')
plot(smanipos_intersect_sos(:,1), smanipos_intersect_sos(:,3),'-m')



%%




