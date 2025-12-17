%{
The model runs the XB model only to obtain the force-calcium curve.
%}

clear all 
addpath XB_model\

exp_Ca = [4.5:.1:5.8,5.85:.05:6.2,6.3:.1:7.0];

%% Baseline 
for i = 1:length(exp_Ca)

% Global Parameters 
data.gpars.ODE_TOL = 1e-8; 

% Metabolite inputs 
data.MgATP_cyto = 8; data.MgADP_cyto = 0; data.Pi_cyto = 0;
% data.MgATP_cyto = 7.8517; data.MgADP_cyto = 0.0514; data.Pi_cyto = 1.2705;

% Sarcomere length inputs 
data.Ls = 2.2; 

% Ca input
data.Ca_i = (10^6)*10^(-exp_Ca(i)); 

% Get parameters and initial conditions 
[init,pars] = parameters_xb(data); 

% Set to tspan to solve for 20 beats
T     = 0.18; dt    = 0.001;  tspan = 0:dt:20*T;
data.T = T; data.dt = dt; data.tspan = tspan; 

% Solve model 
outputs = model_sol(data,pars,init);

sigma_act_RV = outputs.stressess.sigma_act_RV; 

sigma_act_RV_vec(i,:) = sigma_act_RV(end);

% exp_Ca_vec(i,:) = exp_Ca; 

end

F_Nx_max = max(sigma_act_RV_vec);    % Max force, in this case at pCa 4.5 

F_norm = sigma_act_RV_vec./F_Nx_max; % Normalize force by max force
F_50 = F_Nx_max/2;                   % 50% of max force 

% Interpolation to find the pCa value at 50% max force (F_50) i.e. pCa50

% Find the y values of the points above and below F_50 in the force vector
F_top_vec = sigma_act_RV_vec(sigma_act_RV_vec > F_50);     
F_bottom_vec = sigma_act_RV_vec(sigma_act_RV_vec < F_50);
F_1 = F_top_vec(end);
F_2 = F_bottom_vec(1); 

% Get the index to find the corresponding x value in the exp_Ca vector 
x_1 = exp_Ca(sigma_act_RV_vec == F_1);
x_2 = exp_Ca(sigma_act_RV_vec == F_2);

% Interpolate 
exp_Ca_50 = x_1 + (x_2-x_1)/(F_2-F_1)*(F_50-F_1);

% Plot 
figure(100)
clf
hold on 
h1 = plot(-exp_Ca,F_norm,'*-','color',[0.07,.77,0.44]);
h2 = plot(-exp_Ca_50,max(F_norm)/2,'k*');

figure(200)
clf
hold on 
h3 = plot(-exp_Ca,sigma_act_RV_vec,'*-','color',[0.07,.77,0.44]);
h4 = plot(-exp_Ca_50,F_50,'k*');

%% HX 
for i = 1:length(exp_Ca)

% Increase k_on
% data.eta_k_on = 1.5; 

% Decrease k_off 
data.eta_k_off = 0.6;

% Increase kstiff2
data.eta_kstiff = 1.7;

% Increase cooperativity 
% data.eta_K_coop = 2;

% Decrease K_D 
% data.eta_K_D = 0.01; 

% Increase K_T 
% data.eta_K_T = 25; 

% Decrease K_Pi 
% data.eta_K_Pi = .01;

% Increase k_passive 
% data.eta_k_passive = 2; 

% Increase kd
% data.eta_kd = 200; 

% Decrease k3
% data.eta_k3 = 0.35; 

% Increase ka
% data.eta_ka = 2; 

% Increase k1
% data.eta_k1 = 100; 

% Increase km1
% data.eta_km1 = 100; 
% 
% Decrease k2
% data.eta_k2 = 0.1; 
% 
% Increase km2 
% data.eta_km2 = 2.43; 

% Decrease kforce
% data.eta_kforce = .01; 

% Decrease ksr 
% data.eta_ksr = .1; 

% Decrease kmsr 
% data.eta_kmsr = 50; 

% Increase alpha1 
% data.eta_alpha1 = 100; 

% Increase alpha1 
% data.eta_alpha2 = 100; 

% Increase alpha3 
% data.eta_alpha3 = 100; 

% Increase s3
% data.eta_s3 = 100;

% Global Parameters 
data.gpars.ODE_TOL = 1e-8; 

% Metabolite inputs 
data.MgATP_cyto = 8; data.MgADP_cyto = 0; data.Pi_cyto = 0;
% data.MgATP_cyto = 7.8517; data.MgADP_cyto = 0.0514; data.Pi_cyto = 1.2705;
% data.MgATP_cyto = 8*.95; data.MgADP_cyto = 0.05*1.05; data.Pi_cyto = 1.5*1.05;

% Sarcomere length inputs  
data.Ls = 2.2; 

% Ca input
data.Ca_i = (10^6)*10^(-exp_Ca(i)); 

% Get parameters and initial conditions 
[init,pars] = parameters_xb(data); 

% Set to tspan to solve for 20 beats
T     = 0.18; dt    = 0.001;  tspan = 0:dt:20*T;
data.T = T; data.dt = dt; data.tspan = tspan; 

% Solve model 
outputs = model_sol(data,pars,init);

sigma_act_RV = outputs.stressess.sigma_act_RV; 
sigma_act_RV_vec(i,:) = sigma_act_RV(end);

sigma_RV = outputs.stressess.sigma_RV; 
sigma_RV_vec(i,:) = sigma_RV(end);

% exp_Ca_vec(i,:) = exp_Ca; 

end

F_norm = sigma_act_RV_vec./F_Nx_max; % normalize to Nx force

F_Hx_max = max(sigma_act_RV_vec);    % Max force, in this case at pCa 4.5 
F_50 = F_Hx_max/2;                   % 50% of max force 

% Interpolation to find the pCa value at 50% max force (F_50) i.e. pCa50

% Find the y values of the points above and below F_50 in the force vector
F_top_vec = sigma_act_RV_vec(sigma_act_RV_vec > F_50);     
F_bottom_vec = sigma_act_RV_vec(sigma_act_RV_vec < F_50);
F_1 = F_top_vec(end);
F_2 = F_bottom_vec(1); 

% Get the index and find the corresponding y value in the exp_Ca vector 
x_1 = exp_Ca(sigma_act_RV_vec == F_1);
x_2 = exp_Ca(sigma_act_RV_vec == F_2);

% Interpolate 
exp_Ca_50 = x_1 + (x_2-x_1)/(F_2-F_1)*(F_50-F_1);

figure(100)
hold on 
h5 = plot(-exp_Ca,F_norm,'o-','color',[0.6,0.4,1.0]);
h6 = plot(-exp_Ca_50,max(F_norm)/2,'ko');
legend([h1,h5],'Nx','Hx','Location','southeast')
% legend('baseline','kstiff x 1.7','Location','southeast')
xlabel('pCa')
ylabel('F/F_{Nx, max}')
title('RV')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',15)

figure(200)
hold on 
h7 = plot(-exp_Ca,sigma_act_RV_vec,'o-','color',[0.6,0.4,1.0]);
h8 = plot(-exp_Ca_50,F_50,'ko');
legend([h3,h7],{'Nx','Hx'},'Location','southeast')
% legend('baseline','kstiff x 1.7','Location','southeast')
xlabel('pCa')
ylabel('F (kPa)')
title('RV')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',15)
