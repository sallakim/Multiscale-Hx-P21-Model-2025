%{
This script runs the mutliscale cardiovascular model for the selected animals.
%}

close all; clear;clc;

%% Initialization 
% Flags
flag_optpars = 1;      % = 1 will run the optimized parameters, else nominal  
flag_save_figures = 0; % = 1 will save the figures as .eps and .png files
flag_save_outputs = 0; % = 1 will save the outputs as .mat files 

% Select animals to run 
selected_nx = [11 12 51 52 54 55 56];
selected_hx = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62];   

% Animal information 
all_nx_animal_ids = [11 12 51 52 54 55 56];
all_hx_animal_ids = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62];  
nx_labels = ["1", "2", "3", "4", "5", "6", "7"];
hx_labels = ["8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"];

addpath parameters\
addpath opt_pars\
addpath model\
addpath data\

%% Global parameters 
ODE_TOL = 1e-8; 

%% Run energetics model 
TAN_sham = 0.0076;          % 0.0080; % 0.0071–0.0086 mole*(l cell)^(-1)
CRtot_sham = 0.030;         % 0.0267–0.0330 mole*(l cell)^(-1)
TEP_sham = 0.0260;          % 0.0247–0.0298 mole*(l cell)^(-1)
Ox_capacity_sham = 1.0;     % 0.834–1.1526 unitless 
x_ATPase_sham = 1.9/1000;   % 0.9301–1.392 mole sec^(-1) (l cell)^(-1) 

data.TAN_sham = TAN_sham; 
data.CRtot_sham = CRtot_sham; 
data.TEP_sham = TEP_sham; 
data.Ox_capacity_sham = Ox_capacity_sham;
data.x_ATPase_sham = x_ATPase_sham; 

mito  = EnergeticsModelScript(data);

ATP = round(mito(1)); ADP = round(mito(2),2); Pi = round(mito(10),1);

%% Run cardiovascular mechanics  

% Nx 
if exist('selected_nx','var')
    hfig4a = figure(400); clf;
    hfig4c = figure(402); clf;
for n = 1:length(selected_nx)
    % Initialize 
    animal_id = selected_nx(n); % select animals 
    loc = find(animal_id == all_nx_animal_ids);
    nx_label = nx_labels(loc);
    nx_or_hx_flag = 0; % = 0 for nx, = 1 for hx

    % Load data 
    table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
    data = make_datastructure_P21(animal_id,table);
    data.eta_Vtot = 1;
    data.gpars.ODE_TOL = ODE_TOL;      
    data.MgATP_cytoplasm = ATP; 
    data.MgADP_cytoplasm = ADP; 
    data.Pi_cyto         = Pi; 

    % Get parameters 
    optpars_filename = sprintf('opt_pars_Nx%d.mat', animal_id);
    if exist(optpars_filename,'file') == 2 && flag_optpars == 1 
        load(optpars_filename);
        pars(INDMAP) = xopt;
    else
        [pars,~,~,data] = parameters_Nx(data);
    end

    % Run model 
    outputs = model_sol(pars,data);

    % Save outputs 
    if flag_save_outputs == 1
        folder = fullfile(pwd,'outputs');
        filename = sprintf('outputs_Nx%d.mat', animal_id);
        fullpath = fullfile(folder,filename);
        save(fullpath,'outputs','pars','data','animal_id')
    end

    % Plot 
    data_filename = sprintf('Nx%d_data.mat',animal_id);
    load(data_filename);
    % plot_PV_loops(outputs,n,nx_label,nx_or_hx_flag,V_LV_mult,V_LV_mean,P_LV_mult,P_LV_mean,V_RV_mult,V_RV_mean,P_RV_mult,P_RV_mean)
    plot_pv_timecourse(outputs,n,nx_label,nx_or_hx_flag,V_LV_mult,V_LV_mean,P_LV_mult,P_LV_mean,V_RV_mult,V_RV_mean,P_RV_mult,P_RV_mean)
end
end

% Hx 
if exist('selected_hx','var')
    hfig4b = figure(401); clf;
    hfig4d = figure(403); clf;
for n = 1:length(selected_hx)
    % Initialize 
    animal_id = selected_hx(n); % select animals 
    loc = find(animal_id == all_hx_animal_ids);
    hx_label = hx_labels(loc);
    nx_or_hx_flag = 1; % = 0 for nx, = 1 for hx

    % Load data 
    table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
    data = make_datastructure_P21(animal_id,table);
    data.eta_Vtot = 1;
    data.gpars.ODE_TOL = ODE_TOL;      
    data.MgATP_cytoplasm = ATP; 
    data.MgADP_cytoplasm = ADP; 
    data.Pi_cyto         = Pi; 
    data.hx_flag = 1; 

    % Get parameters 
    optpars_filename = sprintf('opt_pars_Hx%d.mat', animal_id);
    if exist(optpars_filename,'file') == 2 && flag_optpars == 1 
        load(optpars_filename);
        pars(INDMAP) = xopt;
    else
        [pars,~,~,data] = parameters_Hx(data);
    end

    % Run model 
    outputs = model_sol(pars,data);

    % Save outputs 
    if flag_save_outputs == 1
        folder = fullfile(pwd,'outputs');
        filename = sprintf('outputs_Hx%d.mat', animal_id);
        fullpath = fullfile(folder,filename);
        save(fullpath,'outputs','pars','data','animal_id')
    end

    % Plot 
    data_filename = sprintf('Hx%d_data.mat',animal_id);
    load(data_filename);
    plot_PV_loops(outputs,n,hx_label,nx_or_hx_flag,V_LV_mult,V_LV_mean,P_LV_mult,P_LV_mean,V_RV_mult,V_RV_mean,P_RV_mult,P_RV_mean)
    plot_pv_timecourse(outputs,n,hx_label,nx_or_hx_flag,V_LV_mult,V_LV_mean,P_LV_mult,P_LV_mean,V_RV_mult,V_RV_mean,P_RV_mult,P_RV_mean)
end
end


%% Save figures 
if flag_save_figures == 1
    if exist('selected_nx','var')
        print(hfig4a,'-dpng',strcat('Figures/','/F4_nx_rv.png'))
        print(hfig4a,'-depsc2',strcat('Figures/','/F4_nx_rv.eps'))
        print(hfig4c,'-dpng',strcat('Figures/','/F4_nx_lv.png'))
        print(hfig4c,'-depsc2',strcat('Figures/','/F4_nx_lv.eps'))
    end
    if exist('selected_hx','var')
        print(hfig4b,'-dpng',strcat('Figures/','/F4_hx_rv.png'))
        print(hfig4b,'-depsc2',strcat('Figures/','/F4_hx_rv.eps'))
        print(hfig4d,'-dpng',strcat('Figures/','/F4_hx_lv.png'))
        print(hfig4d,'-depsc2',strcat('Figures/','/F4_hx_lv.eps'))
    end
end
