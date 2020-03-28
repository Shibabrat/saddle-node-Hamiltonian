%   SCRIPT to compute unstable periodic orbit and its invariant manifolds
%   for the 2 DoF saddle-node Hamiltonian (= TE + PE)
%--------------------------------------------------------------------------
%   Potential energy surface notations:
%
%     Saddle (EQNUM=1)          Well (stable, EQNUM = 2)    
%
%--------------------------------------------------------------------------
% Shibabrat Naik (28-April-2019)
% par = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];

global eqNum deltaE parameters

% Setting up parameters and global variables
N = 4;          % dimension of phase space
MASS_A = 1.0; MASS_B = 1.0; % Mass-weighted momenta

% MU = 0.1;
% ALPHA = 0.05;
% OMEGA = 1.0;
% EPSILON = 1.5;

MU = 4.0;
OMEGA = 3.0;
ALPHA = 1.0;
% EPSILON = 1e-20;
EPSILON = 5.0;

% MU = 0.1;
% ALPHA = 2.0;
% OMEGA = 1.0;
% EPSILON = 1e-20;
% % EPSILON = 0.1;
% EPSILON = 0.125;
% 
% MU = 0.25;
% ALPHA = 2.0;  
% OMEGA = 1.25;
% EPSILON = 0.25;
% 
% MU = 0.25;
% ALPHA = 2.0;
% OMEGA = 1.25;
% EPSILON = 1.0;
% 
% MU = 0.25;
% ALPHA = 2.0;
% OMEGA = 1.25;
% EPSILON = 1.5;

% n_mfd_traj = 1000;
n_mfd_traj = 15;

parameters = [MASS_A MASS_B MU ALPHA OMEGA EPSILON];

eqNum = 1;  
[eqPt] = get_eq_pts_saddlenode2dof(eqNum, parameters);

eSaddle = get_TE_saddlenode2dof(eqPt', parameters); % energy of the saddle eq pt


% [q,x] = SNode_wellLocation(MU, ALPHA, OMEGA, EPSILON)

%% 

nFam = 100; % use nFam = 10 for low energy

% first two amplitudes for continuation procedure to get p.o. family
Ax1  = 2.e-4; % initial amplitude (1 of 2) values to use: 2.e-3
Ax2  = 2*Ax1; % initial amplitude (2 of 2)

tic;

%  get the initial conditions and periods for a family of periodic orbits
po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum),'_saddlenode2dof.txt'];
[po_x0Fam,po_tpFam] = get_POFam_saddlenode2dof(eqNum, Ax1, Ax2, ...
                            nFam, po_fam_file, parameters) ; 

poFamRuntime = toc;

x0podata = [po_x0Fam, po_tpFam];


%%

% begins with a family of periodic orbits and steps until crossing the
% initial condition with target energy 
% fileName = 'x0po_T_energy_case1_L41.txt';
% fileName = 'x0po_T.txt';
deltaE = 0.5;

po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum),'_saddlenode2dof.txt'];
eTarget = eSaddle + deltaE; 
fprintf('Loading the periodic orbit family from data file %s \n',po_fam_file); 
x0podata = importdata(po_fam_file);

po_brac_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                '_brac',num2str(deltaE),'_saddlenode2dof.txt'];
tic;
% [x0poTarget,TTarget] = bracket_POEnergy_bp(eTarget, x0podata, po_brac_file);
[x0poTarget,TTarget] = poBracketEnergy_saddlenode2dof(eTarget, x0podata, ...
                        po_brac_file, parameters);
poTarE_runtime = toc;

save(['model_parameters_eqPt',num2str(eqNum), ...
    '_DelE',num2str(deltaE),'_saddlenode2dof.txt'], ...
    'parameters', '-ASCII', '-double');



%%

% target specific periodic orbit
% Target PO of specific energy with high precision; does not work for the
% model 

po_target_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                    '_DelE',num2str(deltaE),'_saddlenode2dof.txt'];
                
[x0_PO, T_PO, e_PO] = poTargetEnergy_saddlenode2dof(x0poTarget, ...
                        eTarget,po_target_file,parameters);


% data_path = ['./data/UPOs-deltaE',num2str(deltaE),'/x0po_T_energyPO_eqPt', ...
%     num2str(eqNum),'_DelE',num2str(deltaE),'.txt'];
data_path = ['./x0po_T_energyPO_eqPt', num2str(eqNum), ...
    '_DelE', num2str(deltaE), '_saddlenode2dof.txt'];
% data_path = ['./x0po_T_energyPO_eqPt', num2str(eqNum), '_brac', ...
%     num2str(deltaE), '.txt']

x0po = importdata(data_path);

TPOFam = x0po(:,5); 
ePOFam = x0po(:,6);
nMed = size(x0po,1);
tmfd = 5*TPOFam(nMed);
% tmfd = 3.2*TPOFam(nMed);
% tmfd = 6*TPOFam(nMed);

% n_mfd_traj = 25;
frac = 0.1;
del = 1e-6;


%% Stable manifold, negative branch

tmfd = 10.5*TPOFam(nMed);
% deltaE = 0.125;
% stbl = 1; dir = 1;
stbl = -1; dir = -1;

tic;    

[xW,x0W] = get_POManiLocal_saddlenode2dof(x0po(nMed,1:4),TPOFam(nMed),frac, ...
                                            stbl,dir,del,tmfd, ...
                                            n_mfd_traj,parameters);

maniRuntime = toc

energyTube = ePOFam(nMed) ;
% title(['Total energy: ', num2str(energyTube)], ...
%     'interpreter','Latex','Fontsize', 16);
% view(18,9)
view(18,15)
box on
grid on

%% Stable manifold, positive branch

% n_mfd_traj = 1000;
stbl = -1; dir = 1;
% tmfd = 5*TPOFam(nMed);
tmfd = 12.5*TPOFam(nMed);

tic;    

[xW,x0W] = get_POManiLocal_saddlenode2dof(x0po(nMed,1:4),TPOFam(nMed),frac, ...
                                            stbl,dir,del,tmfd, ...
                                            n_mfd_traj,parameters);

maniRuntime = toc

energyTube = ePOFam(nMed) ;
% title(['Total energy: ', num2str(energyTube)], ...
%     'interpreter','Latex','Fontsize', 16);
view(20,25)

%% Unstable manifold, negative branch

stbl = 1; dir = -1;
tmfd = 10.5*TPOFam(nMed);
% tmfd = 2.5*TPOFam(nMed);

tic;    

[xW,x0W] = get_POManiLocal_saddlenode2dof(x0po(nMed,1:4),TPOFam(nMed),frac, ...
                                            stbl,dir,del,tmfd, ...
                                            n_mfd_traj,parameters);

maniRuntime = toc

energyTube = ePOFam(nMed) ;
% title(['Total energy: ', num2str(energyTube)], ...
%     'interpreter','Latex','Fontsize', 16);



%% Unstable manifold, positive branch

% n_mfd_traj = 25;
stbl = 1; dir = 1;
% tmfd = 8*TPOFam(nMed);
tmfd = 11.5*TPOFam(nMed);

tic;    

[xW,x0W] = get_POManiLocal_saddlenode2dof(x0po(nMed,1:4),TPOFam(nMed),frac, ...
                                            stbl,dir,del,tmfd, ...
                                            n_mfd_traj,parameters);

maniRuntime = toc

energyTube = ePOFam(nMed) ;
% title(['Total energy: ', num2str(energyTube)], ...
%     'interpreter','Latex','Fontsize', 16);
view(20,25)


%%


% xx = importdata(['upo_',po_target_file]);
% 
% 
% plot3([xx(:,1);flip(xx(:,1))], [xx(:,2);flip(xx(:,2))], ...
%                 [xx(:,4);flip(-xx(:,4))],'-k','Linewidth',10);
% grid on
% box on
% hold on
% 
% %%
% 
fs = draw_energysurf_saddlenode2dof(deltaE,parameters,0.3);






%%

















