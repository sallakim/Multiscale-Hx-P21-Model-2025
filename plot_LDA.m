clear; clc; 

addpath opt_pars\
addpath model\

nx_animal_ids = [11 12 51 52 54 55 56];
hx_animal_ids = [1 4 5 7 8 9 10 57 58 59 61 62]; 

green = [0.07,.77,0.44]; 
purple = [0.6,0.4,1.0]; 

for i = 1:length(nx_animal_ids)
    animal_id = nx_animal_ids(i);
    filename2 = sprintf('opt_pars_Nx%d.mat', animal_id); % get optimized parameters 
    table = readtable('P21_data_input.xlsx','PreserveVariableNames',true); % get data 
    load(filename2,'xopt')
    opt_pars = exp(xopt); 
    opt_pars_nx(i,:) = opt_pars; 

    loc = find(table.AnimalID == animal_id);
    t_current = table(loc,:);
    ESP_LV_nx(i) = t_current.LVESP;
    EDP_LV_nx(i) = t_current.LVEDP;
    
    ESP_RV_nx(i) = t_current.RVESP;
    EDP_RV_nx(i) = t_current.RVEDP;
    
    ESV_LV_nx(i) = t_current.LVESV;

    ESV_RV_nx(i) = t_current.RVESV;

    SV_LV_nx(i) = t_current.StrokeVolume;

end

for i = 1:length(hx_animal_ids)
    animal_id = hx_animal_ids(i);
    filename2 = sprintf('opt_pars_Hx%d.mat', animal_id); % get optimized parameters 
    table = readtable('P21_data_input.xlsx','PreserveVariableNames',true); % get data 
    load(filename2,'xopt')
    opt_pars = exp(xopt); 
    opt_pars_hx(i,:) = opt_pars; % each row is a separate animal 

    loc = find(table.AnimalID == animal_id);
    t_current = table(loc,:);
    ESP_LV_hx(i) = t_current.LVESP;
    EDP_LV_hx(i) = t_current.LVEDP;
    
    ESP_RV_hx(i) = t_current.RVESP;
    EDP_RV_hx(i) = t_current.RVEDP;
    
    ESV_LV_hx(i) = t_current.LVESV;

    ESV_RV_hx(i) = t_current.RVESV;

    SV_LV_hx(i) = t_current.StrokeVolume;

end


%% Format data for the LDA 
Nx_data = [ESV_LV_nx;ESP_LV_nx;EDP_LV_nx;ESV_RV_nx;ESP_RV_nx;EDP_RV_nx;SV_LV_nx]'; % Rows: animals, Columns: data
Hx_data = [ESV_LV_hx;ESP_LV_hx;EDP_LV_hx;ESV_RV_hx;ESP_RV_hx;EDP_RV_hx;SV_LV_hx]';
Nx_data_length = length(Nx_data); 
Hx_data_length = length(Hx_data); 
Y = [zeros(Nx_data_length,1); ones(Hx_data_length,1)];

X_data = [Nx_data;Hx_data];
X_data_norm = zscore(X_data); % Normalize due to different units / scales 

X_model = [opt_pars_nx;opt_pars_hx];
X_model_norm = zscore(X_model);

X_m_d_norm = [X_data_norm,X_model_norm];

%% Do LDA 
% Do LDA on all parameters (combined data + model)
lda_data = fitcdiscr(X_data_norm, Y);

% Project the high-dimensional data into the LDA space
X_lda_scores_data = X_data_norm * lda_data.Coeffs(1,2).Linear + lda_data.Coeffs(1,2).Const;

% Compute the canonical LDA scores (this gives a multi-D projection)
[~, ldaScores_data, ~] = predict(lda_data, X_data_norm);



% Do LDA on all parameters (combined data + model)
lda_model = fitcdiscr(X_model_norm, Y);

% Project the high-dimensional data into the LDA space
X_lda_scores_model = X_model_norm * lda_model.Coeffs(1,2).Linear + lda_model.Coeffs(1,2).Const;

% Compute the canonical LDA scores (this gives a multi-D projection)
[~, ldaScores_model, ~] = predict(lda_model, X_model_norm);



% Do LDA on all parameters (combined data + model)
lda_m_d = fitcdiscr(X_m_d_norm, Y);

% Project the high-dimensional data into the LDA space
X_lda_scores_m_d = X_m_d_norm * lda_m_d.Coeffs(1,2).Linear + lda_m_d.Coeffs(1,2).Const;

% Compute the canonical LDA scores (this gives a multi-D projection)
[~, ldaScores_m_d, ~] = predict(lda_m_d, X_m_d_norm);




%% Plot
fig100 = figure(100);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

marker_size = 100;               % consistent circle size
jitter_strength = 0.01;         % same jitter amount for all panels

% Data Alone 
nexttile 
hold on

scores_data = X_lda_scores_data;
y_jitter_data = jitter_strength * randn(length(scores_data), 1);

% Plot Nx (green) and Hx (purple)
scatter(scores_data(Y==0), y_jitter_data(Y==0), marker_size, ...
    'filled', 'MarkerFaceColor', green, 'MarkerEdgeColor', 'k');
scatter(scores_data(Y==1), y_jitter_data(Y==1), marker_size, ...
    'filled', 'MarkerFaceColor', purple, 'MarkerEdgeColor', 'k');

% Add numeric labels (1–19) inside circles
for i = 1:length(scores_data)
    text(scores_data(i), y_jitter_data(i), num2str(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', 'w', 'FontWeight', 'bold', 'FontSize', 8);
end

xlabel('LDA score');
title('Data Alone');
ylim([-0.05 0.05]);
yticks([]);
grid on;


% Model Parameters Alone 
nexttile 
hold on

scores_model = X_lda_scores_model;
y_jitter_model = jitter_strength * randn(length(scores_model), 1);

scatter(scores_model(Y==0), y_jitter_model(Y==0), marker_size, ...
    'filled', 'MarkerFaceColor', green, 'MarkerEdgeColor', 'k');
scatter(scores_model(Y==1), y_jitter_model(Y==1), marker_size, ...
    'filled', 'MarkerFaceColor', purple, 'MarkerEdgeColor', 'k');

xlabel('LDA score');
title('Model Parameters Alone');
ylim([-0.05 0.05]);
yticks([]);
grid on;


% Data + Model Parameters 
nexttile 
hold on

scores_m_d = X_lda_scores_m_d;
y_jitter_m_d = jitter_strength * randn(length(scores_m_d), 1);

scatter(scores_m_d(Y==0), y_jitter_m_d(Y==0), marker_size, ...
    'filled', 'MarkerFaceColor', green, 'MarkerEdgeColor', 'k');
scatter(scores_m_d(Y==1), y_jitter_m_d(Y==1), marker_size, ...
    'filled', 'MarkerFaceColor', purple, 'MarkerEdgeColor', 'k');

xlabel('LDA score');
title('Data + Model Parameters');
legend({'Nx','Hx'}, 'Location','best');
ylim([-0.05 0.05]);
yticks([]);
grid on;

%% Plot
figure(1)
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile 
hold on

% Plot Nx subjects (label 0)
scatter(X_lda_scores_data(Y==0), zeros(sum(Y==0),1), 'o', 'filled','MarkerFaceColor',green);

% Plot Hx subjects (label 1)
scatter(X_lda_scores_data(Y==1), zeros(sum(Y==1),1), 'o', 'filled','MarkerFaceColor',purple);

xlabel('LDA score');
title('Data Alone');
ylim([-0.1 0.1]); % Flat 1D line
grid on

nexttile 
hold on

% Plot Nx subjects (label 0)
scatter(X_lda_scores_model(Y==0), zeros(sum(Y==0),1), 'o', 'filled','MarkerFaceColor',green);

% Plot Hx subjects (label 1)
scatter(X_lda_scores_model(Y==1), zeros(sum(Y==1),1), 'o', 'filled','MarkerFaceColor',purple);

xlabel('LDA score');
title('Model Parameters Alone');
ylim([-0.1 0.1]); % Flat 1D line
grid on

nexttile 
hold on

% Plot Nx subjects (label 0)
scatter(X_lda_scores_m_d(Y==0), zeros(sum(Y==0),1), 'o', 'filled','MarkerFaceColor',green);

% Plot Hx subjects (label 1)
scatter(X_lda_scores_m_d(Y==1), zeros(sum(Y==1),1), 'o', 'filled','MarkerFaceColor',purple);

xlabel('LDA score');
title('Data + Model Parameters');
legend('Nx', 'Hx');
ylim([-0.1 0.1]); % Flat 1D line
grid on



%%

print(fig100,'-depsc2',strcat('Figures/','/lda.eps'))
print(fig100,'-dpng',strcat('Figures/','/lda_raw.png'))

