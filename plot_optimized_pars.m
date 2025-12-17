%% Plot optimized parameters boxplots 

clear 
clc

addpath opt_pars\

savefigon = 0; % = 1 will save the figures as .eps 

selected_nx = [11 12 51 52 54 55 56];
selected_hx = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62];  % animal indexes to loop

f_ids = [4 5 9 55 56 57 62]; % indexes of the female animals 

INDMAP = [5, 6, 11, 12, 13, 14, 15, 16, 17, 18]; % optimized paraemter indexes

% params =  {'C_{SA}', 'C_{SV}', 'C_{PA}', 'C_{PV}', ...
%            'R_{SA}', 'R_{PA}', ...
%            'R_m', 'R_a', 'R_t', 'R_p', ...
%            'Amref_{LV}', 'Amref_{SEP}', 'Amref_{RV}', ...
%            'Vw_{LV}', 'Vw_{SEP}', 'Vw_{RV}', ...
%            'k_{TS}', 'k_{TR}'};          

params = {'$C_{SA}$', '$C_{SV}$', '$C_{PA}$', '$C_{PV}$', ...
          '$R_{SA}$', '$R_{PA}$', ...
          '$R_m$', '$R_a$', '$R_t$', '$R_p$', ...
          '$A_{m,ref,LV}$', '$A_{m,ref,SEP}$', '$A_{m,ref,RV}$', ...
          '$V_{w,LV}$', '$V_{w,SEP}$', '$V_{w,RV}$', ...
          '$k_{TS}$','$k_{TR}$'}; % formatted for latex

counter = 0; 

%% Build a matrix with all the optimized parameters 
% Nx 
for i = 1:length(selected_nx)
    animal_id = selected_nx(i);  % get the current animal number 
    
    % Check if the file exists
    filename = sprintf('opt_pars_Nx%d_NEW_largebounds_2025.mat', animal_id);  % get the file name
    if exist(filename, 'file') == 2  
        load(filename);  % load file  
        counter = counter + 1; 
        optpars_mat(:, counter) = exp(xopt);  % save parameters in a matrix 
        % Separate F and M animals based on condition
        if ismember(animal_id, f_ids) % check if the animal is female
            optpars_F_Nx(:, counter) = exp(xopt);
        else
            optpars_M_Nx(:, counter) = exp(xopt);
        end
    end
end
nx_end = counter(end);

% Hx (repeat) 
for i = 1:length(selected_hx)
    animal_id = selected_hx(i);  
    
    % Check if the file exists 
    filename = sprintf('opt_pars_Hx%d_NEW_largebounds_2025.mat', animal_id);
    if exist(filename, 'file') == 2  
        load(filename);  
        counter = counter + 1;  
        optpars_mat(:, counter) = exp(xopt);  
        % Separate F and M animals based on condition
        if ismember(animal_id, f_ids) 
            optpars_F_Hx(:, counter) = exp(xopt);
        else
            optpars_M_Hx(:, counter) = exp(xopt);
        end
    end
end

% Remove empty columns with all zeros
optpars_F_Nx(:,all(optpars_F_Nx == 0))=[];
optpars_F_Hx(:,all(optpars_F_Hx == 0))=[];

optpars_M_Nx(:,all(optpars_M_Nx == 0))=[];
optpars_M_Hx(:,all(optpars_M_Hx == 0))=[];


%% Plots 
x1 = ones(1,length(selected_nx));  % create the x data needed to overlay the swarmchart on the boxchart
x2 = ones(1,length(selected_hx));  % create the x data needed to overlay the swarmchart on the boxchart

% Box plots of the optimized parameters separated by CONDITION 
hfig1 = figure(1);
clf

% Create a tiled layout for 5 row and 2 columns
tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(INDMAP)
    param_index = INDMAP(i); % get the parameter 
    Nx_grp = optpars_mat(i,1:nx_end); % select the corresponding parameters 
    Hx_grp = optpars_mat(i,nx_end+1:end);
    grp = [ones(size(Nx_grp)), 2*ones(size(Hx_grp))]; % Accounts for different Nx, Hx sample sizes 
    plot_opts = optpars_mat(i,:);
    [h,p]  = ttest2(Nx_grp,Hx_grp);
    % p_vec(i) = p; % store the p and h values in a vector 
    % h_vec(i) = h; 
    
    nexttile
    hold on
    boxplot(plot_opts,grp)
    if h == 1
       significance(plot_opts,p)
    end
    % swarmchart(1*x1,plot_opts(1:7),20,[0.07,.77,0.44],'jitter','on','jitterAmount',1e-100,'markerfacecolor',[0.07,.77,0.44])
    % swarmchart(2*x2,plot_opts(8:19),20,[0.6,0.4,1.0],'jitter','on','jitterAmount',1e-100,'markerfacecolor',[0.6,0.4,1.0])

    plotsettings(params,param_index,plot_opts,p) % apply plot settings (below)

    xticklabels({'Nx','Hx'})
end


% Boxplots of the optimized parameters separated by SEX and CONDITION
hfig2 = figure(2);
clf

tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(INDMAP)
    param_index = INDMAP(i);
    grp = [ones(size(optpars_F_Hx(i,:))), 2*ones(size(optpars_M_Hx(i,:)))];
    plot_opts = [optpars_F_Hx(i,:),optpars_M_Hx(i,:)];
    [h,p]  = ttest2(optpars_F_Hx(i,:),optpars_M_Hx(i,:));

    nexttile
    hold on 
    sgtitle('Hx','FontSize',12)
    boxplot(plot_opts,grp)
    if h == 1
        significance(plot_opts,p)
    end

    plotsettings(params,param_index,plot_opts,p)

    xticklabels({'F','M'})
    
end

hfig3 = figure(3);
clf

t = tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(INDMAP)
    param_index = INDMAP(i);
    grp = [ones(size(optpars_F_Nx(i,:))), 2*ones(size(optpars_M_Nx(i,:)))];
    plot_opts = [optpars_F_Nx(i,:),optpars_M_Nx(i,:)];
    % not enough data for t-test

    nexttile
    hold on
    plotsettings(params,param_index,plot_opts,p)
    sgtitle('Nx','FontSize',12)  
    boxplot(plot_opts,grp) 
    xticklabels({'F','M'})
end


% Boxplots of the optimized parameters to compare within the sex between
% groups

hfig4 = figure(4);
clf

tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(INDMAP)
    param_index = INDMAP(i);
    grp = [ones(size(optpars_F_Nx(i,:))), 2*ones(size(optpars_F_Hx(i,:)))];
    plot_opts = [optpars_F_Nx(i,:),optpars_F_Hx(i,:)];
    % not enough data for t-test

    nexttile
    hold on 
    sgtitle('F','FontSize',12)
    boxplot(plot_opts,grp)
    if h == 1
        significance(plot_opts,p)
    end

    plotsettings(params,param_index,plot_opts,p)

    xticklabels({'Nx','Hx'})
    
end

hfig5 = figure(5);
clf

t = tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(INDMAP)
    param_index = INDMAP(i);
    grp = [ones(size(optpars_M_Nx(i,:))), 2*ones(size(optpars_M_Hx(i,:)))];
    plot_opts = [optpars_M_Nx(i,:),optpars_M_Hx(i,:)];
    [h,p]  = ttest2(optpars_M_Nx(i,:),optpars_M_Hx(i,:));
    

    nexttile
    hold on
    plotsettings(params,param_index,plot_opts,p)
    sgtitle('M','FontSize',12)  
    boxplot(plot_opts,grp) 
    xticklabels({'Nx','Hx'})
end


%% Save Figure 
if savefigon == 1
    print(hfig1,'-depsc2',strcat('Figures/','/optimized_parameters.eps'))
end

%% Functions

% Plot significance 
function significance(plot_opts,p)
if 0.01 < p && p<= 0.05
    plot(1.5, max(plot_opts)*0.95, '*k')
elseif 0.001 < p && p<= 0.01
    plot(1.5, max(plot_opts)*0.95, '*k')
    plot(1.5, max(plot_opts)*0.80, '*k')
elseif 0.0001 < p && p<= 0.001
    plot(1.5, max(plot_opts)*0.95, '*k')
    plot(1.5, max(plot_opts)*0.80, '*k')
    plot(1.5, max(plot_opts)*0.65, '*k')
elseif p < 0.0001 
    plot(1.5, max(plot_opts)*0.95, '*k')
    plot(1.5, max(plot_opts)*0.80, '*k')
    plot(1.5, max(plot_opts)*0.65, '*k')
    plot(1.5, max(plot_opts)*0.50, '*k')
end
end

% Reused plot settings 
function plotsettings(params,param_index,plot_opts,p)
set(gcf, 'Renderer', 'painters')
set(gca,'FontSize',12)
title(params(param_index),'interpreter', 'latex','FontSize',20)
xlim([0.5,2.5])
ymin = round(min(plot_opts)*0.8,3);
ymax = round(max(plot_opts)*1.1,2);  
if p<= 0.05
    ymax = round(max(plot_opts)*1.3,2); 
end
ylim([ymin,ymax])
end


%%
% [mean(optpars_mat(1,1:nx_end)),mean(optpars_mat(1,nx_end+1:end))]
% [mean(optpars_mat(2,1:nx_end)),mean(optpars_mat(2,nx_end+1:end))]
% [mean(optpars_mat(3,1:nx_end)),mean(optpars_mat(3,nx_end+1:end))]
% [mean(optpars_mat(4,1:nx_end)),mean(optpars_mat(4,nx_end+1:end))]
% [mean(optpars_mat(5,1:nx_end)),mean(optpars_mat(5,nx_end+1:end))]
% [mean(optpars_mat(6,1:nx_end)),mean(optpars_mat(6,nx_end+1:end))]
% [mean(optpars_mat(7,1:nx_end)),mean(optpars_mat(7,nx_end+1:end))]
% [mean(optpars_mat(8,1:nx_end)),mean(optpars_mat(8,nx_end+1:end))]
% [mean(optpars_mat(9,1:nx_end)),mean(optpars_mat(9,nx_end+1:end))]
% [mean(optpars_mat(10,1:nx_end)),mean(optpars_mat(10,nx_end+1:end))]

% [std(optpars_mat(1,1:nx_end)), std(optpars_mat(1,nx_end+1:end))]
% [std(optpars_mat(2,1:nx_end)), std(optpars_mat(2,nx_end+1:end))]
% [std(optpars_mat(3,1:nx_end)), std(optpars_mat(3,nx_end+1:end))]
% [std(optpars_mat(4,1:nx_end)), std(optpars_mat(4,nx_end+1:end))]
% [std(optpars_mat(5,1:nx_end)), std(optpars_mat(5,nx_end+1:end))]
% [std(optpars_mat(6,1:nx_end)), std(optpars_mat(6,nx_end+1:end))]
% [std(optpars_mat(7,1:nx_end)), std(optpars_mat(7,nx_end+1:end))]
% [std(optpars_mat(8,1:nx_end)), std(optpars_mat(8,nx_end+1:end))]
% [std(optpars_mat(9,1:nx_end)), std(optpars_mat(9,nx_end+1:end))]
% [std(optpars_mat(10,1:nx_end)), std(optpars_mat(10,nx_end+1:end))]

% [h,p] = ttest2(optpars_mat(1,1:nx_end), optpars_mat(1,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(1,1:nx_end),  optpars_mat(1,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(2,1:nx_end),  optpars_mat(2,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(3,1:nx_end),  optpars_mat(3,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(4,1:nx_end),  optpars_mat(4,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(5,1:nx_end),  optpars_mat(5,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(6,1:nx_end),  optpars_mat(6,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(7,1:nx_end),  optpars_mat(7,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(8,1:nx_end),  optpars_mat(8,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(9,1:nx_end),  optpars_mat(9,nx_end+1:end))
% [h,p] = ttest2(optpars_mat(10,1:nx_end), optpars_mat(10,nx_end+1:end))
% 
% [mean(optpars_F_Hx(1,:)),  mean(optpars_M_Hx(1,:))]
% [mean(optpars_F_Hx(2,:)),  mean(optpars_M_Hx(2,:))]
% [mean(optpars_F_Hx(3,:)),  mean(optpars_M_Hx(3,:))]
% [mean(optpars_F_Hx(4,:)),  mean(optpars_M_Hx(4,:))]
% [mean(optpars_F_Hx(5,:)),  mean(optpars_M_Hx(5,:))]
% [mean(optpars_F_Hx(6,:)),  mean(optpars_M_Hx(6,:))]
% [mean(optpars_F_Hx(7,:)),  mean(optpars_M_Hx(7,:))]
% [mean(optpars_F_Hx(8,:)),  mean(optpars_M_Hx(8,:))]
% [mean(optpars_F_Hx(9,:)),  mean(optpars_M_Hx(9,:))]
% [mean(optpars_F_Hx(10,:)), mean(optpars_M_Hx(10,:))]
% 
% [std(optpars_F_Hx(1,:)),  std(optpars_M_Hx(1,:))]
% [std(optpars_F_Hx(2,:)),  std(optpars_M_Hx(2,:))]
% [std(optpars_F_Hx(3,:)),  std(optpars_M_Hx(3,:))]
% [std(optpars_F_Hx(4,:)),  std(optpars_M_Hx(4,:))]
% [std(optpars_F_Hx(5,:)),  std(optpars_M_Hx(5,:))]
% [std(optpars_F_Hx(6,:)),  std(optpars_M_Hx(6,:))]
% [std(optpars_F_Hx(7,:)),  std(optpars_M_Hx(7,:))]
% [std(optpars_F_Hx(8,:)),  std(optpars_M_Hx(8,:))]
% [std(optpars_F_Hx(9,:)),  std(optpars_M_Hx(9,:))]
% [std(optpars_F_Hx(10,:)), std(optpars_M_Hx(10,:))]
% 
% [h,p] = ttest2(optpars_F_Hx(1,:),  optpars_M_Hx(1,:))
% [h,p] = ttest2(optpars_F_Hx(2,:),  optpars_M_Hx(2,:))
% [h,p] = ttest2(optpars_F_Hx(3,:),  optpars_M_Hx(3,:))
% [h,p] = ttest2(optpars_F_Hx(4,:),  optpars_M_Hx(4,:))
% [h,p] = ttest2(optpars_F_Hx(5,:),  optpars_M_Hx(5,:))
% [h,p] = ttest2(optpars_F_Hx(6,:),  optpars_M_Hx(6,:))
% [h,p] = ttest2(optpars_F_Hx(7,:),  optpars_M_Hx(7,:))
% [h,p] = ttest2(optpars_F_Hx(8,:),  optpars_M_Hx(8,:))
% [h,p] = ttest2(optpars_F_Hx(9,:),  optpars_M_Hx(9,:))
% [h,p] = ttest2(optpars_F_Hx(10,:), optpars_M_Hx(10,:))
% 
% [mean(optpars_F_Nx(1,:)),  mean(optpars_M_Nx(1,:))]
% [mean(optpars_F_Nx(2,:)),  mean(optpars_M_Nx(2,:))]
% [mean(optpars_F_Nx(3,:)),  mean(optpars_M_Nx(3,:))]
% [mean(optpars_F_Nx(4,:)),  mean(optpars_M_Nx(4,:))]
% [mean(optpars_F_Nx(5,:)),  mean(optpars_M_Nx(5,:))]
% [mean(optpars_F_Nx(6,:)),  mean(optpars_M_Nx(6,:))]
% [mean(optpars_F_Nx(7,:)),  mean(optpars_M_Nx(7,:))]
% [mean(optpars_F_Nx(8,:)),  mean(optpars_M_Nx(8,:))]
% [mean(optpars_F_Nx(9,:)),  mean(optpars_M_Nx(9,:))]
% [mean(optpars_F_Nx(10,:)), mean(optpars_M_Nx(10,:))]
% 
% [std(optpars_F_Nx(1,:)),  std(optpars_M_Nx(1,:))]
% [std(optpars_F_Nx(2,:)),  std(optpars_M_Nx(2,:))]
% [std(optpars_F_Nx(3,:)),  std(optpars_M_Nx(3,:))]
% [std(optpars_F_Nx(4,:)),  std(optpars_M_Nx(4,:))]
% [std(optpars_F_Nx(5,:)),  std(optpars_M_Nx(5,:))]
% [std(optpars_F_Nx(6,:)),  std(optpars_M_Nx(6,:))]
% [std(optpars_F_Nx(7,:)),  std(optpars_M_Nx(7,:))]
% [std(optpars_F_Nx(8,:)),  std(optpars_M_Nx(8,:))]
% [std(optpars_F_Nx(9,:)),  std(optpars_M_Nx(9,:))]
% [std(optpars_F_Nx(10,:)), std(optpars_M_Nx(10,:))]
% 
% [h,p] = ttest2(optpars_M_Nx(1,:),  optpars_M_Hx(1,:))
% [h,p] = ttest2(optpars_M_Nx(2,:),  optpars_M_Hx(2,:))
% [h,p] = ttest2(optpars_M_Nx(3,:),  optpars_M_Hx(3,:))
% [h,p] = ttest2(optpars_M_Nx(4,:),  optpars_M_Hx(4,:))
% [h,p] = ttest2(optpars_M_Nx(5,:),  optpars_M_Hx(5,:))
% [h,p] = ttest2(optpars_M_Nx(6,:),  optpars_M_Hx(6,:))
% [h,p] = ttest2(optpars_M_Nx(7,:),  optpars_M_Hx(7,:))
% [h,p] = ttest2(optpars_M_Nx(8,:),  optpars_M_Hx(8,:))
% [h,p] = ttest2(optpars_M_Nx(9,:),  optpars_M_Hx(9,:))
% [h,p] = ttest2(optpars_M_Nx(10,:), optpars_M_Hx(10,:))
