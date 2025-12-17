%{
This plots the sensitivity analysis results. 
%}
%% Initialize 
clear; clc;
addpath sens_results\

savefigs = 0; % = 1 saves figures

% Define Nx and Hx subject IDs
nx_animal_ids = [11 12 51 52 54 55 56];
hx_animal_ids = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62];

subjectIDs = [nx_animal_ids,hx_animal_ids];

params = {'$C_{SA}$', '$C_{SV}$', '$C_{PA}$', '$C_{PV}$', ...
          '$R_{SA}$', '$R_{PA}$', ...
          '$R_m$', '$R_a$', '$R_t$', '$R_p$', ...
          '$A_{m,ref,LV}$', '$A_{m,ref,SEP}$', '$A_{m,ref,RV}$', ...
          '$V_{w,LV}$', '$V_{w,SEP}$', '$V_{w,RV}$', ...
          '$k_{TS}$','$k_{TR}$'}; % latex format

params =  {'C_{SA}', 'C_{SV}', 'C_{PA}', 'C_{PV}', ...
           'R_{SA}', 'R_{PA}', ...
           'R_m', 'R_a', 'R_t', 'R_p', ...
           'A_{m,ref,LV}', 'A_{m,ref,SEP}', 'A_{m,ref,RV}', ...
           'V_{w,LV}', 'V_{w,SEP}', 'V_{w,RV}', ...
           'k_{TS}', 'k_{TR}'};     

color_codes = [0.29,0.49,0.27;
    0.1,0.8,0.3;
    0.7,1.0,0.0;
    0.0,1.0,1.0;   
    0.0,0.6,1.0;
    0.1,0.0,1.0;
    0.6,0.0,1.0;
    1.0,0.0,1.0;
    1.0,0.0,0.0;
    1.0,0.4,0.0;
    0.71,0.54,0.13;
    0.6,0.1,0.2;
    ];

%% Prepare the data (functions defined at the bottom)
% Process animals and load data 
process_animals(nx_animal_ids, 'Nx');
process_animals(hx_animal_ids, 'Hx');

% Loop through each subject and compute the scalar
for i = 1:length(subjectIDs)
    % Get the variable name for the current subject
    varName = sprintf('sens%d', subjectIDs(i));
    
    % Check if the variable exists in the workspace 
    if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
        scalarName = sprintf('sens_scalar_%d', subjectIDs(i)); % create variable name
        evalin('base', sprintf('%s = sqrt(sum(%s.^2));', scalarName, varName)); % computer
    else
        warning('Variable %s does not exist in the workspace.', varName);
    end
end

% Create matrices with the scalar sensitivities for Nx and Hx animals
sens_matrix_nx = create_sens_matrix(nx_animal_ids);
sens_matrix_hx = create_sens_matrix(hx_animal_ids);

% Rank scalar sensitivities and reorder the parameters 
[~, sens_order_hx] = sort(mean(sens_matrix_hx),'descend');
params_sorted_hx = params(sens_order_hx);
sens_matrix_hx = sens_matrix_hx(:,sens_order_hx);

[~, sens_order_nx] = sort(mean(sens_matrix_nx),'descend');
params_sorted_nx = params(sens_order_nx);
sens_matrix_nx = sens_matrix_nx(:,sens_order_nx);

%% Plotting 
% Assign plot colors 
green = [0.07,.77,0.44];
purple = [0.6,0.4,1.0];

% Create plot settings 
function applyPlotSettings(params_sorted)
    ylim([0, 1]);
    yticks(0:.2:1);
    xticks(1:length(params_sorted));
    xticklabels(params_sorted);
    xaxisproperties = get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'tex'; 
    xtickangle(90)
    pbaspect([1 1 1]);
    

    set(gcf, 'Renderer', 'painters');
    set(gca, 'FontSize', 12,'FontName','Cambria');
end

% Create x data for swarmchart 
x_nx = (1:length(params)) .* ones(size(sens_matrix_nx, 1), 1);
x_hx = (1:length(params)) .* ones(size(sens_matrix_hx, 1), 1);

hfig2 = figure(2);
clf
set(hfig2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.7]);  % Larger figure


% First subplot (Nx)
ax1=subplot(1,2,1);
hold on
title('Nx','FontSize',20)
yline(.1, '--', 'LineWidth', .5); % Dashed line at y = 5
boxplot(sens_matrix_nx)
swarmchart(x_nx, sens_matrix_nx, 25, green,'filled')%,'MarkerFaceAlpha',0.5)
applyPlotSettings(params_sorted_nx)
set(ax1, 'Position', [0.07 0.15 0.45 0.75]);  % [left bottom width height]

% Second subplot (Hx)
ax2=subplot(1,2,2);
hold on
title('Hx','FontSize',20)
yline(.1, '--', 'LineWidth', .5); % Dashed line at y = 5
boxplot(sens_matrix_hx)
swarmchart(x_hx, sens_matrix_hx, 25, purple,'filled')%,'MarkerFaceAlpha',0.5)
applyPlotSettings(params_sorted_hx)
set(ax2, 'Position', [0.55 0.15 0.45 0.75]);  % Shift right, same height/width


%% Printing 
if savefigs == 1
    print(hfig2,'-depsc2',strcat('Figures/','/F2.eps'))
end

%% Functions Used 

% Normalize each section of the 'sens' matrix by dividing by the corresponding mean data 'y0'
function sens = normalize_sens(sens, y0)
    indices = [1, 330, 659, 988, length(sens) + 1];  % Define the section breakpoints
    for i = 1:4
        sens(indices(i):indices(i+1)-1, :) = sens(indices(i):indices(i+1)-1, :) / mean(y0(indices(i):indices(i+1)-1));
    end
end

% Function to process animals
function process_animals(animal_ids, group)
    for i = 1:length(animal_ids)
        % Construct the filename dynamically
        fileName = sprintf('sens_Rat_%s%d.mat', group, animal_ids(i));
        
        % Load the data if the file exists
        if exist(fileName, 'file')
            load(fileName);
            
            % Normalize the 'sens' matrix
            sens = normalize_sens(sens, y0);
            
            % Generate variable name (e.g., 'sens11')
            varName = sprintf('sens%d', animal_ids(i));
            
            % Store the normalized 'sens' into the dynamically named variable
            assignin('base', varName, sens);
        else
            warning('File %s does not exist', fileName);
        end
    end
end

% Function to create a sens matrix
function sens_matrix = create_sens_matrix(animal_ids)
    sens_matrix = [];
    for i = 1:length(animal_ids)
        varName = sprintf('sens_scalar_%d', animal_ids(i));
        sens_data = evalin('base', varName);
        sens_matrix = [sens_matrix; sens_data / max(sens_data)];
    end
end