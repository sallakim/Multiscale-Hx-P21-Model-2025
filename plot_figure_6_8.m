% This script plots cardiac power (Figure 6) and correlations (Figure 8)
clear 
clc

addpath opt_pars\
addpath outputs\

savefigs = 1; % = 1 saves the figures

selected_nx = [11 12 51 52 54 55 56];
selected_hx = [1 4 5 7 8 9 10 57 58 59 61 62]; 

%% Get Data
for i = 1:length(selected_nx)
    nx_or_hx_flag = 0; % = 0 for Nx, = 1 for Hx
    animal_id = selected_nx(i); % get the current animal 
    filename = sprintf('outputs_Nx%d.mat', animal_id); % get the files to load
    filename2 = sprintf('opt_pars_Nx%d.mat', animal_id);
    load(filename)
    load(filename2,'pars_opt')
    pars_opt = exp(pars_opt);
    
    [time, ...
        sigma_XB_LV, sigma_XB_SEP, sigma_XB_RV, ...
        Lsc_LV, Lsc_SEP, Lsc_RV, ...
        P_RV, P_PA, ...
        V_LV, V_RV, ...
        Cm_LV, Cm_SEP, Cm_RV, ...
        Q_a_valve, Q_m_valve] = load_outputs(animal_id,nx_or_hx_flag);
    
    Q_a_valve = Q_a_valve(1:floor(length(Q_a_valve)/2)); 
    Q_m_valve = Q_m_valve(1:floor(length(Q_m_valve)/2)); 

    work_LV = cumtrapz(Lsc_LV,sigma_XB_LV)*10^(-6);
    work_SEP = cumtrapz(Lsc_SEP,sigma_XB_SEP)*10^(-6);
    work_RV = cumtrapz(Lsc_RV,sigma_XB_RV)*10^(-6);

    power_LV_nx(i) = max(abs(diff(work_LV)./diff(time)));
    power_SEP_nx(i) = max(abs(diff(work_SEP)./diff(time)));
    power_RV_nx(i) = max(abs(diff(work_RV)./diff(time)));
    
    SW_RV_nx(i) = trapz(P_RV,V_RV);
    
    max_strain_SEP_nx (i) = max((Lsc_SEP - Lsc_SEP(1))./Lsc_SEP(1)); 
    max_strain_LV_nx (i) = max((Lsc_LV - Lsc_LV(1))./Lsc_LV(1)); 
    
    R_PA_nx(i) = pars_opt(6);
    
    P_RV_max_nx(i) = max(P_RV);
    
    P_PA_mean_nx(i) = mean(P_PA); 
    
    Amref_RV_nx(i) = pars_opt(13);

    C_PA_nx(i) = pars_opt(3);
    
    % Find ejection
    i_AVO = find(diff(Q_a_valve) > 0,1,'first'); % Start of sytole 
    i_AVC = find(diff(Q_a_valve) < 0,1,'last');  % End of systole 
    
    % Find filling
    i_MVO = find(diff(Q_m_valve) > 0,1,'first'); % Start of sytole 
    i_MVC = find(diff(Q_m_valve) < 0,1,'last');   % End of diastole 
    
    Cm_systole = Cm_SEP(i_AVO:i_AVC); 
    Cm_systole_nx(i) = mean(Cm_systole); 
    
    Cm_drop_nx(i) = Cm_SEP(1)-min(Cm_SEP(1:i_AVC));

    P_PV = outputs.pressures.P_PV; 
    P_PV_mean_nx(i) = mean(P_PV); 

    L_s_RV_mean_Nx(i) = mean(Lsc_RV);

end

for i = 1:length(selected_hx)
    nx_or_hx_flag = 1; % = 0 for Nx, = 1 for Hx
    animal_id = selected_hx(i); % get the current animal
    filename = sprintf('outputs_Hx%d.mat', animal_id); % get the files to load
    filename2 = sprintf('opt_pars_Hx%d.mat', animal_id);
    load(filename)
    load(filename2,'pars_opt')
    pars_opt = exp(pars_opt);
    
    [time, ...
        sigma_XB_LV, sigma_XB_SEP, sigma_XB_RV, ...
        Lsc_LV, Lsc_SEP, Lsc_RV, ...
        P_RV, P_PA, ...
        V_LV, V_RV, ...
        Cm_LV, Cm_SEP, Cm_RV, ...
        Q_a_valve, Q_m_valve] = load_outputs(animal_id,nx_or_hx_flag);
    
    Q_a_valve = Q_a_valve(1:floor(length(Q_a_valve)/2)); 
    Q_m_valve = Q_m_valve(1:floor(length(Q_m_valve)/2)); 
    
    work_LV = cumtrapz(Lsc_LV,sigma_XB_LV)*10^(-6); % kPa*m = N*m / m^2 = J / m^2
    work_SEP = cumtrapz(Lsc_SEP,sigma_XB_SEP)*10^(-6);
    work_RV = cumtrapz(Lsc_RV,sigma_XB_RV)*10^(-6);

    power_LV_hx(i) = max(abs(diff(work_LV)./diff(time))); % J/s / m^2 = W / m^2 
    power_SEP_hx(i) = max(abs(diff(work_SEP)./diff(time)));
    power_RV_hx(i) = max(abs(diff(work_RV)./diff(time)));

    SW_RV_hx(i) = trapz(P_RV,V_RV);
    
    max_strain_SEP_hx (i) = max((Lsc_SEP - Lsc_SEP(1))./Lsc_SEP(1)); 
    max_strain_LV_hx (i) = max((Lsc_LV - Lsc_LV(1))./Lsc_LV(1)); 
    
    R_PA_hx(i) = pars_opt(6);
    
    P_RV_max_hx(i) = max(P_RV);
    
    P_PA_mean_hx(i) = mean(P_PA); 
    
    Amref_RV_hx(i) = pars_opt(13);

    C_PA_hx(i) = pars_opt(3);
    
    % Find ejection
    i_AVO = find(diff(Q_a_valve) > 0,1,'first'); % Start of sytole 
    i_AVC = find(diff(Q_a_valve) < 0,1,'last');  % End of systole 
    
    % Find filling
    i_MVO = find(diff(Q_m_valve) > 0,1,'first'); % Start of sytole 
    i_MVC = find(diff(Q_m_valve) < 0,1,'last');   % End of diastole 
    
    Cm_systole = Cm_SEP(i_AVO:i_AVC); 
    Cm_systole_hx(i) = mean(Cm_systole); 
    
    Cm_drop_hx(i) = Cm_SEP(1)-min(Cm_SEP(1:i_AVC));

    P_PV = outputs.pressures.P_PV; 
    P_PV_mean_hx(i) = mean(P_PV); 

    L_s_RV_mean_Hx(i) = mean(Lsc_RV);

end

%% Myofiber Power 
power_LV = [power_LV_nx,power_LV_hx];
power_RV = [power_RV_nx,power_RV_hx];

grp = [ones(size(power_RV_nx)), 2*ones(size(power_RV_hx))]; % in order to plot box plots with different amounts of data

x1 = ones(1,length(selected_nx));  % create the x data needed to overlay the swarmchart on the boxchart
x2 = ones(1,length(selected_hx)); % create the x data needed to overlay the swarmchart on the boxchart

hfig1 = figure(1);
clf 

subplot(1,2,1) 
hold on 
boxplot(power_RV,grp,'Colors','k')
swarmchart(1*x1,power_RV(1:7),50,[0.07,.77,0.44],'jitter','on','jitterAmount',0.05,'MarkerFaceColor',[0.07,.77,0.44])
swarmchart(2*x2,power_RV(8:19),50,[0.6,0.4,1.0],'jitter','on','jitterAmount',0.05,'MarkerFaceColor',[0.6,0.4,1.0])
xlim([0.5 2.5])
ylim([20 90])
xticklabels({'Nx','Hx'})
ylabel('Myofiber Power Intensity (W per m^2)')
set(gcf, 'Renderer', 'painters')
set(gca,'fontsize',12)
% set(gca,'TickLabelInterpreter','latex')
title('RV')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
[h,p]  = ttest2(power_RV(1:7),power_RV(8:19));
if h == 1
    significance(xt,yt,p)
end

subplot(1,2,2)
hold on 
boxplot(power_LV,grp,'Colors','k')
swarmchart(1*x1,power_LV(1:7),50,[0.07,.77,0.44],'jitter','on','jitterAmount',0.05,'MarkerFaceColor',[0.07,.77,0.44])
swarmchart(2*x2,power_LV(8:19),50,[0.6,0.4,1.0],'jitter','on','jitterAmount',0.05,'MarkerFaceColor',[0.6,0.4,1.0])
xlim([0.5 2.5])
ylim([10 130])
xticklabels({'Nx','Hx'})
set(gcf, 'Renderer', 'painters')
set(gca,'fontsize',12)
% set(gca,'TickLabelInterpreter','latex')
title('LV')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
[h,p]  = ttest2(power_LV(1:7),power_LV(8:19));
if h == 1
    significance(xt,yt,p)
end

%% Scatter Plot 
function [RHO_S,RHO_P] = corrplot(input1_nx,input1_hx,input2_nx,input2_hx)
green = [0.07,.77,0.44]; 
purple = [0.6,0.4,1.0]; 
input1 = [input1_nx,input1_hx];
input2 = [input2_nx,input2_hx];
[RHO_S,PVAL_S] = corr(input1',input2','Type','Spearman');
[RHO_P,PVAL_P] = corr(input1',input2','Type','Pearson');
[rho_s,pval_s] = corr(input1([1:7,9:10,14:end])',input2([1:7,9:10,14:end])','Type','Spearman');
[rho_p,pval_p] = corr(input1([1:7,9:10,14:end])',input2([1:7,9:10,14:end])','Type','Pearson');
for i = 1:length(input1_nx)
    plot(input1_nx(i),input2_nx(i),'o','Color',green,'MarkerFaceColor',green)
    text(input1_nx(i),input2_nx(i),num2str(i),FontSize=10)
end
for i = 1:length(input2_hx)
    plot(input1_hx(i),input2_hx(i),'o','Color',purple,'MarkerFaceColor',purple)
    text(input1_hx(i),input2_hx(i),num2str(i+length(input1_nx)),FontSize=10)
end
% title({["Spearman, \rho = " + num2str(round(RHO_S,2))], ["Pearson, \rho = " + num2str(round(RHO_P,2))]})
title(["Pearson, \rho = " + num2str(round(RHO_P,2))])

disp(pval_p)
set(gca, 'FontSize', 12);
end

hfig2 = figure(2);
clf

% Create a tiled layout for 1 row and 3 columns
t = tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot (Cm_drop vs Mean PA Pressure) 
nexttile
hold on
corrplot(P_PA_mean_nx,P_PA_mean_hx,Cm_drop_nx,Cm_drop_hx);
ylabel('${\Delta}C_{m,SEP}$','interpreter','latex','FontSize',20)
xlabel('mPAP','FontSize',20)

% Subplot (Cm_drop vs RV CP) 
nexttile
hold on 
corrplot(Cm_drop_nx,Cm_drop_hx,SW_RV_nx,SW_RV_hx);
ylabel('RV SW','FontSize',20)
xlabel('${\Delta}C_{m,SEP}$','interpreter','latex','FontSize',20)

% % Subplot (Amref_RV vs Mean PA Pressure)
% nexttile
% hold on 
% corrplot(P_PA_mean_nx,P_PA_mean_hx,Amref_RV_nx,Amref_RV_hx);
% ylabel('$A_{m,ref,RV}$','interpreter','latex','FontSize',20)
% xlabel('mPAP','FontSize',20)


nexttile
hold on
corrplot(Cm_drop_nx,Cm_drop_hx,Amref_RV_nx,Amref_RV_hx); 
ylabel('$A_{m,ref,RV}$','interpreter','latex','FontSize',20)
xlabel('${\Delta}C_{m,SEP}$','interpreter','latex','FontSize',20)


%% t-tests
[h,p] = ttest2(L_s_RV_mean_Nx,L_s_RV_mean_Hx); % h = 1 reject the hypothesis 
fprintf('mean RV Ls nx vs mean RV Ls hx: \n')
fprintf('nx: %.2d, hx %.2d \n',mean(L_s_RV_mean_Nx),mean(L_s_RV_mean_Hx))
fprintf('std nx: %.2d, std hx %.2d \n',std(L_s_RV_mean_Nx),std(L_s_RV_mean_Hx))
if h == 1
    fprintf('The means are different, p = %.2d \n\n',p)
else
    fprintf('The means are not different, p = %.2d \n\n',p)
end

[h,p] = ttest2(Cm_systole_nx,Cm_systole_hx); % h = 1 reject the hypothesis 
fprintf('Average curvature over systole for nx vs hx: \n')
fprintf('nx: %.2d, hx %.2d \n',mean(Cm_systole_nx),mean(Cm_systole_hx))
fprintf('std nx: %.2d, std hx %.2d \n',std(Cm_systole_nx),std(Cm_systole_hx))
if h == 1
    fprintf('The means are different, p = %.2d \n\n',p)
else
    fprintf('The means are not different, p = %.2d \n\n',p)
end

[h,p] = ttest2(Cm_drop_nx,Cm_drop_hx);
fprintf('Curvature drop during systole for nx vs hx: \n')
fprintf('mean deltaCm_nx = %.2f, mean deltaCm_hx = %.2f\n', mean(Cm_drop_nx),mean(Cm_drop_hx))
fprintf('std deltaCm_nx = %.2f, std deltaCm_hx = %.2f\n', std(Cm_drop_nx),std(Cm_drop_hx))
if h == 1
    fprintf('The means are different, p = %.2d \n\n',p)
else
    fprintf('The means are not different, p = %.2d \n\n',p)
end

[h,p] = ttest2(P_PA_mean_nx,P_PA_mean_hx); % h = 1 reject the hypothesis 
fprintf('mPAP nx vs mPAP hx: \n')
fprintf('nx: %.2d, hx %.2d \n',mean(P_PA_mean_nx),mean(P_PA_mean_hx))
fprintf('std nx: %.2d, std hx %.2d \n',std(P_PA_mean_nx),std(P_PA_mean_hx))
if h == 1
    fprintf('The means are different, p = %.2d \n\n',p)
else
    fprintf('The means are not different, p = %.2d \n\n',p)
end


[h,p] = ttest2(power_RV_nx,power_RV_hx); % h = 1 reject the hypothesis 
fprintf('RV power nx vs RV power hx: \n')
fprintf('nx: %.2d, hx %.2d \n',mean(power_RV_nx),mean(power_RV_hx))
fprintf('std nx: %.2d, std hx %.2d \n',std(power_RV_nx),std(power_RV_hx))
if h == 1
    fprintf('The means are different, p = %.2d \n',p)
else
    fprintf('The means are not different, p = %.2d \n\n',p)
end

% Pulmonary artery compliance 

% Pulmonary Artery Resistance 

% 
R_C_hx = R_PA_hx.*C_PA_hx;
R_C_nx = R_PA_nx.*C_PA_nx;
[h,p] = ttest2(R_C_nx,R_C_hx);

%% Print Figures
if savefigs == 1
    print(hfig1,'-depsc2',strcat('Figures/','/myofiber_power.eps'))
    print(hfig2,'-depsc2',strcat('Figures/','/scatter_plot.eps'))
end

%% Significance function 
function significance(xt,yt,p)
if 0.01 < p && p<= 0.05
    plot(mean(xt([1 2])), max(yt)*0.95, '*k')
elseif 0.001 < p && p<= 0.01
    plot(mean(xt([1 2])), max(yt)*0.95, '*k')
    plot(mean(xt([1 2])), max(yt)*0.85, '*k')
elseif 0.0001 < p && p<= 0.001
    plot(mean(xt([1 2])), max(yt)*0.95, '*k')
    plot(mean(xt([1 2])), max(yt)*0.85, '*k')
    plot(mean(xt([1 2])), max(yt)*0.75, '*k')
elseif p < 0.0001 
    plot(mean(xt([1 2])), max(yt)*0.95, '*k')
    plot(mean(xt([1 2])), max(yt)*0.85, '*k')
    plot(mean(xt([1 2])), max(yt)*0.75, '*k')
    plot(mean(xt([1 2])), max(yt)*0.65, '*k')
end
end