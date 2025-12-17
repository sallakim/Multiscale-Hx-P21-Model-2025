%{
This script runs the local sensitivity analysis for the multiscale 
cardiovascular model to determine the influence of the "pars" in minimizing 
y (the cost function i.e. the difference between the model output and the 
actual value)
%}

clear; 

addpath model\
addpath parameters\
addpath sensitivity\

% n_cores = 4; 
% delete(gcp('nocreate'))
% p = parpool('local',n_cores); 

%% Initialize 
selected_animals = [51 52];
all_nx_animal_ids = [11 12 51 52 54 55 56];
all_hx_animal_ids = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62]; 

params =  {'$C_{SA}$' '$C_{SV}$' '$C_{PA}$'  '$C_{PV}$' ...       
    '$R_{SA}$' '$R_{PA}$' ...
    '$R_m$' '$R_a$' '$R_t$' '$R_p$'   ... 
    '$Amref_{LV}$' '$Amref_{SEP}$' '$Amref_{RV}$' ...
    '$Vw_{LV}$' '$Vw_{SEP}$' '$Vw_{RV}$' ...
    '$k_{TS}$','$k_{TR}$'};                                 % Reduced parameters names for sensitivity analysis 

% params =  {'$C_{SA}$' '$C_{SV}$' '$C_{PA}$'  '$C_{PV}$' ...       
%     '$R_{SA}$' '$R_{PA}$' ...
%     '$R_m$' '$R_a$' '$R_t$' '$R_p$'   ... 
%     '$A_{m,ref,LV}$' '$A_{m,ref,SEP}$' '$A_{m,ref,RV}$' ...
%     '$Vw_{LV}$' '$Vw_{SEP}$' '$Vw_{RV}$' ...
%     '$k_{TS}$','$k_{TR}$' ...
%     '$k_a$', '$k_d$' ...                          
%     '$K_{Pi}$', '$k_1$', '$k_{-1}$', '$\alpha_1$' ...
%     '$k_2$', '$k_{-2}$', '$\alpha_2$'...          
%     '$K_D$', '$K_T$', '$k_3$', '$\alpha_3$', '$s_3$'...
%     '$k_{on}$', '$k_{off,LV}$', '$k_{off,RV}$', '$K_{coop}$' ...
%     '$k_{sr}$', '$k_{msr}$', '$k_{force}$' ...
%     '$k_{passive,LV}$', '$k_{passive,RV}$', '$\gamma$' ...
%     '$k_{stiff1,LV}$', '$k_{stiff1,RV}$', ...
%     '$k_{stiff2,LV}$', '$k_{stiff2,RV}$' ...
%     '$\delta_R$', '$\eta$', '$K_{se}$' ...
%     '$L_{s,ref}$', '$L_{sc0}$', '$L_{thick}$','$L_{hbare}$','$L_{thin}$' ...
%     };   

for n = 1:length(selected_animals)

    animal_id = selected_animals(n);
    
    %% Get Data
    if ismember(animal_id,all_nx_animal_ids) == 1 % Check if the animal is Nx
        table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
        data = make_datastructure_P21(animal_id,table);
    elseif ismember(animal_id,all_hx_animal_ids) == 1 % Chec, if the animal is Hx
        table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
        data = make_datastructure_P21(animal_id,table);
        data.hx_flag = 1;
    end
    
    ODE_TOL  = 1e-8; 
    data.gpars.ODE_TOL = ODE_TOL; 
    data.MgATP_cytoplasm = 8;
    data.MgADP_cytoplasm = 0.05;
    data.Pi_cyto         = 1.3;
    data.eta_Vtot = 1;
    
    %% Get Nominal Parameters 
    if ismember(animal_id,all_nx_animal_ids) == 1
        [pars,~,~,data] = parameters_Nx(data); 
    elseif ismember(animal_id,all_hx_animal_ids) == 1
        [pars,~,~,data] = parameters_Hx(data); 
    end
    
    %% Run Sensitivity Analysis 
    y0 = LSA_wrap(pars,data);
    
    step_size = .01;
    
    [sens,F] = local_sensitivity(pars,@myfun,y0,step_size,data);
    
    %% Plot
    plot_sens(y0,pars,sens,params,animal_id)
    
    %% Save in a .mat
    if ismember(animal_id,all_nx_animal_ids) == 1
        filename = sprintf('sens_Rat_Nx%d.mat', animal_id);
    elseif ismember(animal_id,all_hx_animal_ids) == 1
        filename = sprintf('sens_Rat_Hx%d.mat', animal_id);
    end
    folder = fullfile(pwd,'sens_results');
    fullpath = fullfile(folder,filename);
    save(fullpath,'sens','y0','params','pars','animal_id')

end

% delete(gcp('nocreate'))

%% Built-in function
function y = myfun(pars,data) 
y = LSA_wrap(pars,data);
end 