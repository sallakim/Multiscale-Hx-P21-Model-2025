%{
This script plots the septal curvature, shortening, and sarcomere length
for LV, SEP, and RV.
%}

addpath outputs 

savefigs = 1; % = 1

% Select animals to run 
selected_nx = [11 12 51 52 54 55 56];
selected_hx = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62];  

% Animal information 
all_nx_animal_ids = [11 12 51 52 54 55 56];
all_hx_animal_ids = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62];  
nx_labels = ["1", "2", "3", "4", "5", "6", "7"];
hx_labels = ["8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"];

%% Name Colors 
color_codes = [
    0.29,0.49,0.27;
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

gray   = [.1 .1 .1];
blue   = [0.0,1.0,1.0];   
orange = [1.0,0.4,0.0];

green  = [0.07,.77,0.44]; 
purple = [0.6,0.4,1.0]; 


%% Initialize 
counter = 0; 

fig1 = figure(700); clf;
fig2 = figure(701); clf;
fig3 = figure(702); clf;

y1_patch = [-.2, -.2, -.18, -.18];
y2_patch = [0 0 .25 .25];
y3_patch = [-5 -5 -4.75 -4.75];

%% Nx
if exist('selected_nx','var')

% Initialize matrices 
Cm_SEP_400_Nx = zeros(length(selected_nx), 400);
Ls_LV_400_Nx  = zeros(length(selected_nx), 400);
Ls_SEP_400_Nx = zeros(length(selected_nx), 400);
Ls_RV_400_Nx  = zeros(length(selected_nx), 400);

for i = 1:length(selected_nx)
    counter = counter + i;
    animal_id = selected_nx(i);
    filename = sprintf('outputs_Nx%d.mat', animal_id); % get the filename
    load(filename,'outputs')

    % Load the outputs 
    time = outputs.time;
    
    Ls_LV  = outputs.lengths.Lsc_LV;
    Ls_SEP = outputs.lengths.Lsc_SEP;
    Ls_RV  = outputs.lengths.Lsc_RV;
    
    Cm_LV  = outputs.curvatures.Cm_LV;
    Cm_SEP = outputs.curvatures.Cm_SEP;
    Cm_RV  = outputs.curvatures.Cm_RV;

    T_norm = 2; 
    time_norm = linspace(0,T_norm,length(Ls_LV));
    T = 1; % period 

    Q_a_valve = outputs.flows.Q_a_valve; 
    Q_a_valve = Q_a_valve(1:floor(length(Q_a_valve)/2)); 
    Q_m_valve = outputs.flows.Q_m_valve;
    Q_m_valve = Q_m_valve(1:floor(length(Q_m_valve)/2));
    
    [T_AVO,T_AVC,T_MVC,T_MVO] = getedes(Q_a_valve,Q_m_valve,time_norm);

    % SHORTENING
    % LV
    fig1 = figure(700);
    subplot(2,3,1)
    hold on
    title('LV')
    plot(time_norm,(Ls_LV - Ls_LV(1))./Ls_LV(1),'Color', green) %color_codes(i,:))
    xlim([0,T])
    ylim([-.2,.15])
    set(gca,'FontSize',10)
    makepatch(y1_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % SEP
    subplot(2,3,2)
    hold on
    title('SEP')
    plot(time_norm,(Ls_SEP - Ls_SEP(1))./Ls_SEP(1),'Color', green) %color_codes(i,:))
    xlim([0,T])
    ylim([-.2,.15])
    set(gca,'FontSize',10)
    makepatch(y1_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % RV
    subplot(2,3,3)
    hold on
    title('RV')
    plot(time_norm,(Ls_RV - Ls_RV(1))./Ls_RV(1),'Color', green) %color_codes(i,:))
    xlim([0,T])
    ylim([-.2,.15])
    set(gca,'FontSize',10)
    makepatch(y1_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % CURVATURE 
    % SEP
    fig2 = figure(701);
    subplot(2,3,1)
    hold on 
    title('LV')
    plot(time_norm,Cm_LV,'Color',green) %color_codes(i,:))
    xlim([0,T])
    ylim([-5 -3])
    set(gca,'FontSize',10)
    makepatch(y3_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    subplot(2,3,2)
    hold on 
    title('SEP')
    plot(time_norm,Cm_SEP,'Color',green) %color_codes(i,:))
    xlim([0,T])
    ylim([-.25 4.75])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    subplot(2,3,3)
    hold on 
    title('RV')
    plot(time_norm,Cm_RV,'Color',green) %color_codes(i,:))
    xlim([0,T])
    ylim([-.25 4.75])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % LENGTH
    % LV
    fig3 = figure(702);
    sgtitle('Sarcomere Length') 
    subplot(2,3,1)
    hold on
    title('LV')
    plot(time_norm,(Ls_LV/10000),'Color',green)
    xlim([0,T])
    ylim([1 3])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % SEP
    subplot(2,3,2)
    hold on
    title('Septum')
    plot(time_norm,(Ls_SEP/10000),'Color',green)
    xlim([0,T])
    ylim([1 3])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % RV
    subplot(2,3,3)
    hold on
    title('RV')
    plot(time_norm,(Ls_RV/10000),'Color',green)
    xlim([0,T])
    ylim([1 3])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)
    
    % Normalize by time 
    T_400 = linspace(0,1,400); 

    Ls_LV_400  = interp1(time_norm,Ls_LV/10000,T_400);
    Ls_SEP_400 = interp1(time_norm,Ls_SEP/10000,T_400);
    Ls_RV_400  = interp1(time_norm,Ls_RV/10000,T_400);

    Cm_LV_400  = interp1(time_norm,Cm_LV,T_400); 
    Cm_SEP_400 = interp1(time_norm,Cm_SEP,T_400); 
    Cm_RV_400  = interp1(time_norm,Cm_RV,T_400); 

    eps_LV  = (Ls_LV - Ls_LV(1))./Ls_LV(1); 
    eps_SEP = (Ls_SEP - Ls_SEP(1))./Ls_SEP(1); 
    eps_RV  = (Ls_RV - Ls_RV(1))./Ls_RV(1); 

    eps_LV_400  = interp1(time_norm,eps_LV,T_400); 
    eps_SEP_400 = interp1(time_norm,eps_SEP,T_400); 
    eps_RV_400  = interp1(time_norm,eps_RV,T_400); 

    % Store in a matrix  
    Ls_LV_400_Nx(i,:)   = Ls_LV_400;
    Ls_SEP_400_Nx(i,:)  = Ls_SEP_400;
    Ls_RV_400_Nx(i,:)   = Ls_RV_400;
    Cm_LV_400_Nx(i,:)   = Cm_LV_400; 
    Cm_SEP_400_Nx(i,:)  = Cm_SEP_400; 
    Cm_RV_400_Nx(i,:)   = Cm_RV_400; 
    eps_LV_400_Nx(i,:)  = eps_LV_400; 
    eps_SEP_400_Nx(i,:) = eps_SEP_400; 
    eps_RV_400_Nx(i,:)  = eps_RV_400; 
end
end

%% Hx
if exist('selected_hx','var')

% Initialize
Cm_SEP_400_Hx = zeros(length(selected_hx), 400);
Ls_LV_400_Hx  = zeros(length(selected_hx), 400);
Ls_SEP_400_Hx = zeros(length(selected_hx), 400);
Ls_RV_400_Hx  = zeros(length(selected_hx), 400);

for i = 1:length(selected_hx)
    counter = counter + i; 
    animal_id = selected_hx(i);
    filename = sprintf('outputs_Hx%d.mat', animal_id);
    load(filename,'outputs')

    % Load the outputs 
    time = outputs.time;

    Ls_LV  = outputs.lengths.Lsc_LV;
    Ls_SEP = outputs.lengths.Lsc_SEP;
    Ls_RV  = outputs.lengths.Lsc_RV;
    
    Cm_LV  = outputs.curvatures.Cm_LV;
    Cm_SEP = outputs.curvatures.Cm_SEP;
    Cm_RV  = outputs.curvatures.Cm_RV;

    T_norm = 2;
    time_norm = linspace(0,T_norm,length(Ls_LV));
    T = 1; 
    
    Q_a_valve = outputs.flows.Q_a_valve; 
    Q_a_valve = Q_a_valve(1:floor(length(Q_a_valve)/2)); 
    Q_m_valve = outputs.flows.Q_m_valve;
    Q_m_valve = Q_m_valve(1:floor(length(Q_m_valve)/2));
    
    [T_AVO,T_AVC,T_MVC,T_MVO] = getedes(Q_a_valve,Q_m_valve,time_norm);

    % SHORTENING 
    figure(700);
    subplot(2,3,4)
    hold on
    plot(time_norm,(Ls_LV - Ls_LV(1))./Ls_LV(1),'Color', purple) %color_codes(i,:))
    xlim([0,T])
    ylim([-.2,.15])
    set(gca,'FontSize',10)
    makepatch(y1_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    subplot(2,3,5)
    hold on
    plot(time_norm,(Ls_SEP - Ls_SEP(1))./Ls_SEP(1),'Color', purple) %color_codes(i,:))
    xlim([0,T])
    ylim([-.2,.15])
    set(gca,'FontSize',10)
    makepatch(y1_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    subplot(2,3,6)
    hold on
    plot(time_norm,(Ls_RV - Ls_RV(1))./Ls_RV(1),'Color', purple) %color_codes(i,:))
    xlim([0,T])
    ylim([-.2,.15])
    set(gca,'FontSize',10)
    makepatch(y1_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % CURVATURE
    figure(701);
    subplot(2,3,4)
    hold on 
    plot(time_norm,Cm_LV,'Color', purple) %color_codes(i,:))
    xlim([0,T])
    ylim([-5 -3])
    set(gca,'FontSize',10)
    makepatch(y3_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    subplot(2,3,5)
    hold on 
    plot(time_norm,Cm_SEP,'Color', purple) %color_codes(i,:))
    xlim([0,T])
    ylim([-.25 4.75])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    subplot(2,3,6)
    hold on 
    plot(time_norm,Cm_RV,'Color', purple) %color_codes(i,:))
    xlim([0,T])
    ylim([-.25 4.75])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % LENGTH
    % LV
    figure(702);
    subplot(2,3,4)
    hold on
    title('LV')
    plot(time_norm,(Ls_LV/10000),'Color',purple)
    xlim([0,T])
    ylim([1 3])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % SEP
    subplot(2,3,5)
    hold on
    title('Septum')
    plot(time_norm,(Ls_SEP/10000),'Color',purple)
    xlim([0,T])
    ylim([1 3])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    % RV
    subplot(2,3,6)
    hold on
    title('RV')
    plot(time_norm,(Ls_RV/10000),'Color',purple)
    xlim([0,T])
    ylim([1 3])
    set(gca,'FontSize',10)
    makepatch(y2_patch,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)

    T_400 = linspace(0,1,400); 

    Ls_LV_400  = interp1(time_norm,Ls_LV/10000,T_400);
    Ls_SEP_400 = interp1(time_norm,Ls_SEP/10000,T_400);
    Ls_RV_400  = interp1(time_norm,Ls_RV/10000,T_400);

    Cm_LV_400  = interp1(time_norm,Cm_LV,T_400); 
    Cm_SEP_400 = interp1(time_norm,Cm_SEP,T_400); 
    Cm_RV_400  = interp1(time_norm,Cm_RV,T_400); 

    eps_LV  = (Ls_LV - Ls_LV(1))./Ls_LV(1); 
    eps_SEP = (Ls_SEP - Ls_SEP(1))./Ls_SEP(1); 
    eps_RV  = (Ls_RV - Ls_RV(1))./Ls_RV(1); 

    eps_LV_400  = interp1(time_norm,eps_LV,T_400); 
    eps_SEP_400 = interp1(time_norm,eps_SEP,T_400); 
    eps_RV_400  = interp1(time_norm,eps_RV,T_400); 

    Ls_LV_400_Hx(i,:)   = Ls_LV_400;
    Ls_SEP_400_Hx(i,:)  = Ls_SEP_400;
    Ls_RV_400_Hx(i,:)   = Ls_RV_400;
    Cm_LV_400_Hx(i,:)   = Cm_LV_400; 
    Cm_SEP_400_Hx(i,:)  = Cm_SEP_400; 
    Cm_RV_400_Hx(i,:)   = Cm_RV_400; 
    eps_LV_400_Hx(i,:)  = eps_LV_400; 
    eps_SEP_400_Hx(i,:) = eps_SEP_400; 
    eps_RV_400_Hx(i,:)  = eps_RV_400; 
end
end

%% Plotting the mean and add axis labels
% Nx
if exist('selected_nx','var')
    % ESP Mean
    figure(700)
    subplot(2,3,1); hold on; plot(T_400, mean(eps_LV_400_Nx), 'k', 'LineWidth', 2);
    subplot(2,3,2); hold on; plot(T_400, mean(eps_SEP_400_Nx), 'k', 'LineWidth', 2);
    subplot(2,3,3); hold on; plot(T_400, mean(eps_RV_400_Nx), 'k', 'LineWidth', 2);

    % Curvature Mean
    figure(701)
    subplot(2,3,1); hold on; plot(T_400, mean(Cm_LV_400_Nx), 'k', 'LineWidth', 2);
    subplot(2,3,2); hold on; plot(T_400, mean(Cm_SEP_400_Nx), 'k', 'LineWidth', 2);
    subplot(2,3,3); hold on; plot(T_400, mean(Cm_RV_400_Nx), 'k', 'LineWidth', 2);

    % Length Mean
    figure(702)
    subplot(2,3,1); hold on; plot(T_400, mean(Ls_LV_400_Nx), 'k', 'LineWidth', 2);
    subplot(2,3,2); hold on; plot(T_400, mean(Ls_SEP_400_Nx), 'k', 'LineWidth', 2);
    subplot(2,3,3); hold on; plot(T_400, mean(Ls_RV_400_Nx), 'k', 'LineWidth', 2);
end

% Hx
if exist('selected_hx','var')
    % ESP Mean
    figure(700)
    subplot(2,3,4); hold on; plot(T_400, mean(eps_LV_400_Hx), 'k', 'LineWidth', 2);
    subplot(2,3,5); hold on; plot(T_400, mean(eps_SEP_400_Hx), 'k', 'LineWidth', 2);
    subplot(2,3,6); hold on; plot(T_400, mean(eps_RV_400_Hx), 'k', 'LineWidth', 2);

    % Curvature Mean
    figure(701)
    subplot(2,3,4); hold on; plot(T_400, mean(Cm_LV_400_Hx), 'k', 'LineWidth', 2);
    subplot(2,3,5); hold on; plot(T_400, mean(Cm_SEP_400_Hx), 'k', 'LineWidth', 2);
    subplot(2,3,6); hold on; plot(T_400, mean(Cm_RV_400_Hx), 'k', 'LineWidth', 2);

    % Length Mean
    figure(702)
    subplot(2,3,4); hold on; plot(T_400, mean(Ls_LV_400_Hx), 'k', 'LineWidth', 2);
    subplot(2,3,5); hold on; plot(T_400, mean(Ls_SEP_400_Hx), 'k', 'LineWidth', 2);
    subplot(2,3,6); hold on; plot(T_400, mean(Ls_RV_400_Hx), 'k', 'LineWidth', 2);
end

% Axis labeling 
figure(700)
set(gcf, 'Renderer', 'painters')
han1 = axes(gcf, 'Visible', 'off');
han1.XLabel.Visible = 'on';
han1.YLabel.Visible = 'on';
xlabel(han1, 'Fraction of Cardiac Cycle', 'FontSize', 12);
ylabel(han1, 'Shortening', 'FontSize', 12);

figure(701)
set(gcf, 'Renderer', 'painters')
han2 = axes(gcf, 'Visible', 'off');
han2.XLabel.Visible = 'on';
han2.YLabel.Visible = 'on';
xlabel(han2, 'Fraction of Cardiac Cycle', 'FontSize', 12);
ylabel(han2, 'Curvature (cm^{-1})', 'FontSize', 12);

figure(702)
set(gcf, 'Renderer', 'painters')
han3 = axes(gcf, 'Visible', 'off');
han3.XLabel.Visible = 'on';
han3.YLabel.Visible = 'on';
xlabel(han3, 'Fraction of Cardiac Cycle', 'FontSize', 12);
ylabel(han3, 'Length (\mum)', 'FontSize', 12);

%% Save Figures 
if savefigs == 1
print(fig1,'-dpng',strcat('Figures/','/F6b.png'))
print(fig2,'-dpng',strcat('Figures/','/F6a.png'))

print(fig1,'-depsc2',strcat('Figures/','/F6b.eps'))
print(fig2,'-depsc2',strcat('Figures/','/F6a.eps'))
end

%% Repeating Functions

function makepatch(y,T_AVO, T_AVC, T_MVO, T_MVC, orange, blue)
    x = [T_AVO T_AVC T_AVC T_AVO];
    patch(x,y,orange,'linestyle','none')
    alpha(.1)

    x = [T_MVO T_MVC T_MVC T_MVO];
    patch(x,y,blue,'linestyle','none')
    alpha(.1)
end

function [T_AVO,T_AVC,T_MVC,T_MVO] = getedes(Q_a_valve,Q_m_valve,time_norm)
    % Find ejection
    i_AVO = find(diff(Q_a_valve) > 0,1,'first'); % Start of sytole 
    i_AVC = find(diff(Q_a_valve) < 0,1,'last');  % End of systole 
    
    T_AVO = time_norm(i_AVO); 
    T_AVC = time_norm(i_AVC); 
    
    % Find filling
    i_MVO = find(diff(Q_m_valve) > 0,1,'first'); % Start of diastole 
    i_MVC = find(diff(Q_m_valve) < 0,1,'last');   % End of diastole 
    if i_MVC < i_MVO
        i_MVO = length(i_MVO);
    end
    if i_MVO < i_AVC
        i_MVO = find(diff(Q_m_valve) == 0,1,'last')+1; 
        i_MVC = length(Q_m_valve);
    end
    T_MVO = time_norm(i_MVO); 
    T_MVC = time_norm(i_MVC); 
end

