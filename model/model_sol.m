function [outputs,rout,J] = model_sol(adjpars,data)

%{ 
    This function solves the time-varying in vivo version of the model to
    steady-state and then calculates 2 steady-state beats. 
    Inputs: 
    adjpars         - vector of adjustable parameters 
    data            - input data structure with data and global parameters
    Outputs: 
    outputs         - structure with all pertinent model outputs to plot 
    rout            - residual vector 
    J               - cost functional 
%} 

%% Initialization  

plot_switch = 0; % 0 = off, 1 = on 

% "undo" log from parameters.m
adjpars = exp(adjpars); 

tspan = data.tspan;  
dt    = data.dt; 

T = data.T; 

ODE_TOL = data.gpars.ODE_TOL;

fixpars = data.fixpars;

%% Parameters 

% Unstressed volumes 
V_SA_u = fixpars(10);
V_SV_u = fixpars(11); 
V_PA_u = fixpars(12); 
V_PV_u = fixpars(13); 

% Compliance 
C_SA = adjpars(1); 

%% Get initial conditions

if isfield(data,'eta_Vtot')
    % Field exists 
    eta_Vtot = data.eta_Vtot; 
else
    % Field does not exist 
    eta_Vtot = 1; 
end 

init = initialconditions(adjpars,data,eta_Vtot);

%% Set mass matrix M for DAE 
M = speye(length(init));
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 
% m = 15; % Number of states 
% dd = ones(m,1); 
% dd(1:4) = 0; 
% M = spdiags(dd,0,m,m); 
% opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);

% try

%% Solve model 

% Use a while loop to allow model to converge to steady-state 
ndone = 0; 
while ndone == 0 

    opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);
    sol  = ode15s(@model,[tspan(1) tspan(end)],init,opts,adjpars,data);
    
    % ndone = 1; %%%% can comment out if you want to run ss code

    % Plot intermediate states if plot_switch is "on"
    if plot_switch == 1 
        % Displacements states 1-4
        figure(111)
        clf
        plot(sol.x,sol.y(1:4,:))
        legend('1','2','3','4')
        
        % Sarcomere lengths states 5-7
        figure(112)
        clf
        plot(sol.x,sol.y(5:7,:))
        legend('5','6','7')
        
        % Volumes states LV, SA
        figure(113)
        clf
        plot(sol.x,sol.y(8:9,:))
        legend('8','9')
        
        % Volumes states 
        figure(114)
        clf
        plot(sol.x,sol.y(43:46,:))
        legend('10','11','12','13')
    end 

    if sol.x(end) ~= tspan(end) 
        % Check to see if the model solved to the end of tspan. If not, set the
        % initial conditions for the next loop at the start of the previous
        % beat and solve for 10 beats 
        t = sol.x(1):dt:sol.x(end); 
        beats = mod(t,T); 
        x = find(round(beats,3) == 0);
        y = find(t(x) <= t(end)); 
        tspan = tspan(1):dt:tspan(x(y(end)));
    
        sols = deval(sol,tspan);
        init  = sols(:,x(y(end-1))); 
    
        tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T; 
         
    else 
        % If the model has successfully solved at least 10 beats, then we can
        % assess whether the model has reached steady state 
        
        % Extract sarcomere lengths and systemic arterial pressure (Psa)     
        sols = deval(sol,tspan);
            
        Lsc_LV  = sols(5,:) ;
        Lsc_SEP = sols(6,:) ; 
        Lsc_RV  = sols(7,:) ; 
        P_SA    = sols(9,:) / C_SA * 7.5; 
        
        % Plot extracted quantities if plot_switch is "on" 
        if plot_switch == 1
            figure(101)
            clf
            hold on 
            plot(tspan,Lsc_LV) 
            
            figure(102)
            clf
            hold on 
            plot(tspan,Lsc_SEP) 
            
            figure(103)
            clf
            hold on 
            plot(tspan,Lsc_RV) 
            
            figure(104)
            clf % new fig with each time series
            hold on 
            plot(tspan,P_SA) 
        end
        
        % Find the last 5 beats of the simulation 
        xx = find(tspan >= tspan(end) - 5*T); 
        
        % Set a peak threshold as half of the amplitude of the last 5 beats 
        threshold_LV  = (max(Lsc_LV(xx))  - min(Lsc_LV(xx)))/2;
        threshold_SEP = (max(Lsc_SEP(xx)) - min(Lsc_SEP(xx)))/2;
        threshold_RV  = (max(Lsc_RV(xx))  - min(Lsc_RV(xx)))/2;
        
        % Determine the length of half of the cardiac cycle 
        half_per = round((T/2)/dt); 
        
        % Find peaks for the last 5 beats 
        [pks_Lsc_LV, loc_pks_Lsc_LV]  = findpeaks(...
            Lsc_LV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_LV); 
        [pks_Lsc_SEP,loc_pks_Lsc_SEP] = findpeaks(...
            Lsc_SEP,'MinPeakDistance',half_per,'MinPeakProminence',threshold_SEP); 
        [pks_Lsc_RV, loc_pks_Lsc_RV]  = findpeaks(...
            Lsc_RV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_RV); 
        [pks_P_SA,   loc_pks_P_SA]    = findpeaks(...
            P_SA,'MinPeakDistance',half_per); 
        
        % Exclude the last peak (so there are 4 peaks)
        pks_Lsc_LV  = pks_Lsc_LV(end-5:end-1); 
        pks_Lsc_SEP = pks_Lsc_SEP(end-5:end-1); 
        pks_Lsc_RV  = pks_Lsc_RV(end-5:end-1);
        pks_P_SA    = pks_P_SA(end-5:end-1); 
        
        % Find the locations of the peaks 
        loc_pks_Lsc_LV  = loc_pks_Lsc_LV(end-5:end-1); 
        loc_pks_Lsc_SEP = loc_pks_Lsc_SEP(end-5:end-1); 
        loc_pks_Lsc_RV  = loc_pks_Lsc_RV(end-5:end-1); 
        loc_pks_P_SA    = loc_pks_P_SA(end-5:end-1); 
        
        % Find the times where the peaks occur 
        t_pks_Lsc_LV  = tspan(loc_pks_Lsc_LV);
        t_pks_Lsc_SEP = tspan(loc_pks_Lsc_SEP);
        t_pks_Lsc_RV  = tspan(loc_pks_Lsc_RV);
        t_pks_P_SA    = tspan(loc_pks_P_SA); 
        
        % Create a linear regression through the peaks 
        pf_Lsc_LV  = polyfit(t_pks_Lsc_LV,pks_Lsc_LV,1); 
        pf_Lsc_SEP = polyfit(t_pks_Lsc_SEP,pks_Lsc_SEP,1); 
        pf_Lsc_RV  = polyfit(t_pks_Lsc_RV,pks_Lsc_RV,1); 
        pf_P_SA    = polyfit(t_pks_P_SA,pks_P_SA,1); 
        
        % Extract the slope of the regression line 
        slope_Lsc_LV  = pf_Lsc_LV(1);
        slope_Lsc_SEP = pf_Lsc_SEP(1);
        slope_Lsc_RV  = pf_Lsc_RV(1);
        slope_P_SA    = pf_P_SA(1);
        
        % Plot regression line through peaks if plot_switch is "on" 
        if plot_switch == 1
            % Draw the regression line through the peaks
            y_Lsc_LV  = polyval(pf_Lsc_LV,tspan); 
            y_Lsc_SEP = polyval(pf_Lsc_SEP,tspan); 
            y_Lsc_RV  = polyval(pf_Lsc_RV,tspan); 
            y_P_SA    = polyval(pf_P_SA,tspan); 
            
            % LV sarcomere length 
            figure(101)
            hold on 
            plot(t_pks_Lsc_LV,pks_Lsc_LV,'r*')
            plot(tspan,y_Lsc_LV,'k')
            
            % SEP sarcomere length 
            figure(102)
            hold on 
            plot(t_pks_Lsc_SEP,pks_Lsc_SEP,'r*')
            plot(tspan,y_Lsc_SEP,'k')
            
            % RV sarcomere length 
            figure(103)
            hold on 
            plot(t_pks_Lsc_RV,pks_Lsc_RV,'r*')
            plot(tspan,y_Lsc_RV,'k')
            
            % Systemic arterial pressure 
            figure(104)
            hold on 
            plot(t_pks_P_SA,pks_P_SA,'r*')
            plot(tspan,y_P_SA,'k') 
        end
        
        % If the slope is sufficiently small (i.e. flat), we have reached 
        % steady state 
        slope_lim = 1e-2;
            % Stopping criteria 
            if abs(slope_P_SA) < slope_lim && abs(slope_Lsc_LV) < slope_lim && ...
                    abs(slope_Lsc_SEP) < slope_lim && abs(slope_Lsc_RV) < slope_lim
                ndone = 1; 
            end 
            
        % If we have not reached steady-state, solve the model for 10 more
        % beats and reassess convergence 
        beats = mod(tspan,T); 
        x = find(round(beats,3) == 0);
        y = find(tspan(x) <= tspan(end));
        tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T;
        init  = sols(:,x(y(end-1))); 
    end 
end

% After determining that the model is in steady-state, solve 2 more beats 
time = [0:dt:2*T]; 
sol  = ode15s(@model,[time(1) time(end)],init,opts,adjpars,data);
sols = deval(sol,time);  % sols length (# initial conditions, # of timesteps) 
sols = sols'; 

%% Calculate other time-varying model quantities (pressures, flows, etc.) 

o = zeros(41,length(time));  
for i = 1:length(time) 
    [~,o(:,i)] = model(time(i),sols(i,:),adjpars,data);
end 

%% Create output structure  

outputs.time = time; 

% Convert cm 
displacements.xm_LV  = sols(:,1) ; 
displacements.xm_SEP = sols(:,2) ; 
displacements.xm_RV  = sols(:,3) ; 
displacements.ym     = sols(:,4) ; 

% Convert cm to um
lengths.Lsc_LV  = sols(:,5) * 1e4;
lengths.Lsc_SEP = sols(:,6) * 1e4; 
lengths.Lsc_RV  = sols(:,7) * 1e4; 

% cm^3 = mL
volumes.V_LV = sols(:,8); 
volumes.V_RV = sols(:,9) ; 
volumes.V_SV = (sols(:,10) + V_SV_u); 
volumes.V_PV = (sols(:,11) + V_PV_u); 
volumes.V_SA = (sols(:,12) + V_SA_u); 
volumes.V_PA = (sols(:,13) + V_PA_u); 
volumes.Vtot = (sum(sols(end,8:13)) + V_SA_u + V_SV_u + V_PA_u + V_PV_u); 

% unitless 
probabilities.P_A1_0_LV = sols(:,14);
probabilities.P_A2_0_LV = sols(:,17);
probabilities.P_A3_0_LV = sols(:,20);
probabilities.N_LV = sols(:,23);
probabilities.P0_LV = o(39,:);
probabilities.P_A1_0_SEP = sols(:,25);
probabilities.P_A2_0_SEP = sols(:,28);
probabilities.P_A3_0_SEP = sols(:,31);
probabilities.N_SEP = sols(40,:);
probabilities.P0_SEP = o(:,33);
probabilities.P_A1_0_RV = sols(:,36);
probabilities.P_A2_0_RV = sols(:,39);
probabilities.P_A3_0_RV = sols(:,42);
probabilities.N_RV = sols(:,45);
probabilities.P0_RV = o(41,:);

% Convert kPa to mmHg
pressures.P_LV = o(1,:) * 7.5; 
pressures.P_SA = o(2,:) * 7.5; 
pressures.P_SV = o(3,:) * 7.5; 
pressures.P_RV = o(4,:) * 7.5; 
pressures.P_PA = o(5,:) * 7.5; 
pressures.P_PV = o(6,:) * 7.5; 

% cm^3
wallvolumes.Vm_LV  = o(7,:)  ; 
wallvolumes.Vm_SEP = o(8,:) ; 
wallvolumes.Vm_RV  = o(9,:) ; 

% cm^2
areas.Am_LV  = o(10,:) ; 
areas.Am_SEP = o(11,:) ; 
areas.Am_RV  = o(12,:) ; 

% cm^(-1)
curvatures.Cm_LV  = o(13,:) ;
curvatures.Cm_SEP = o(14,:) ;
curvatures.Cm_RV  = o(15,:) ; 

strains.eps_LV  = o(16,:); 
strains.eps_SEP = o(17,:); 
strains.eps_RV  = o(18,:); 

stresses.passive.sigma_pas_LV  = o(19,:);
stresses.passive.sigma_pas_SEP = o(20,:);
stresses.passive.sigma_pas_RV  = o(21,:);

stresses.active.sigma_XB_LV  = o(22,:);
stresses.active.sigma_XB_SEP = o(23,:);
stresses.active.sigma_XB_RV  = o(24,:);

stresses.total.sigma_LV  = o(25,:);
stresses.total.sigma_SEP = o(26,:);
stresses.total.sigma_RV  = o(27,:);

% Convert cm^3 s^(-1) = mL min^(-1)
flows.Q_m_valve = o(28,:) * 60; 
flows.Q_a_valve = o(29,:) * 60; 
flows.Q_t_valve = o(30,:) * 60; 
flows.Q_p_valve = o(31,:) * 60; 

flows.Q_SA = o(32,:) * 60; 
flows.Q_PA = o(33,:) * 60; 

tensions.Tm_LV  = o(34,:);
tensions.Tm_SEP = o(35,:); 
tensions.Tm_RV  = o(36,:); 

pressures.P_peri = o(37,:) * 7.5'; 

outputs.XB.Ca_i = o(38,:);

outputs.volumes       = volumes; 
outputs.pressures     = pressures; 
outputs.displacements = displacements; 
outputs.areas         = areas;
outputs.wallvolumes   = wallvolumes; 
outputs.curvatures    = curvatures; 
outputs.strains       = strains; 
outputs.stresses      = stresses;
outputs.lengths       = lengths; 
outputs.flows         = flows; 
outputs.tensions      = tensions; 
outputs.probabilities   = probabilities;

%% Sensitivity Analysis 

% rout = []; 

%% Chamber Volumes 


if isfield(data,'V_LV_avg') == 1 
    V_LV_avg = data.V_LV_avg;
    V_RV_avg = data.V_RV_avg;
    P_LV_avg = data.P_LV_avg;
    P_RV_avg = data.P_RV_avg;
else
    V_LV_avg = 1; 
    V_RV_avg = 1;
    P_LV_avg = 1;
    P_RV_avg = 1;
end


ESV_LV_data = data.ESV_LV;
EDV_LV_data = data.EDV_LV;
ESV_RV_data = data.ESV_RV;
EDV_RV_data = data.EDV_RV;

ESP_LV_data = data.ESP_LV*7.5;
EDP_LV_data = data.EDP_LV*7.5;
ESP_RV_data = data.ESP_RV*7.5;
EDP_RV_data = data.EDP_RV*7.5;

% discritize time to match number elements in pressure/volume data
T_data = linspace(0,T,50);

% make model output vectors same length as data 
 
sol  = ode15s(@model,[T_data(1) T_data(end)],init,opts,adjpars,data);
sols = deval(sol,T_data);  % sols length (# initial conditions, # of timesteps) 
sols = sols'; 

o = zeros(41,length(T_data));  
for i = 1:length(T_data) 
    [~,o(:,i)] = model(T_data(i),sols(i,:),adjpars,data);
end 

V_LV_model = sols(:,8); 
P_LV_model = o(1,:) * 7.5; 
V_RV_model = sols(:,9) ; 
P_RV_model = o(4,:) * 7.5; 

[EDV_model, EDP_model, ESV_model, ESP_model] = getEDESvals(V_LV_model, V_RV_model, P_LV_model, P_RV_model);


rout_1 = (V_LV_avg/1000 - V_LV_model)./(max(V_LV_avg)/1000);
rout_2 = (V_RV_avg/1000 - V_RV_model)./(max(V_RV_avg)/1000);
rout_3 = (P_LV_avg - P_LV_model')./(max(P_LV_avg));
rout_4 = (P_RV_avg - P_RV_model')./(max(P_RV_avg));

rout_5 = (EDV_LV_data - EDV_model(1))./EDV_LV_data; 
rout_6 = (EDV_RV_data - EDV_model(2))./EDV_RV_data; 

rout_7 = (ESV_LV_data - ESV_model(1))./ESV_LV_data; 
rout_8 = (ESV_RV_data - ESV_model(2))./ESV_RV_data; 

rout_9  = (EDP_LV_data - EDP_model(1))./EDP_LV_data; 
rout_11 = (EDP_RV_data - EDP_model(2))./EDP_RV_data; 

rout_10 = (ESP_LV_data - ESP_model(1))./ESP_LV_data; 
rout_12 = (ESP_RV_data - ESP_model(2))./ESP_RV_data; 


rout = [1/sqrt(50)*rout_1; 1/sqrt(50)*rout_2; 1/sqrt(50)*rout_3; 1/sqrt(50)*rout_4;
    rout_5; rout_6; rout_7; rout_8;
    rout_9; rout_10; rout_11; rout_12];

% figure(1000)
% subplot(2,1,1)
% hold on 
% plot(V_LV_model, P_LV_model)
% plot(V_LV_avg/1000,P_LV_avg,'k', 'linewidth',2)
% subplot(2,1,2)
% hold on 
% plot(V_RV_model, P_RV_model)
% plot(V_RV_avg/1000,P_RV_avg,'k', 'linewidth',2)
% set(gca,'FontSize',20)


% 
% catch 
%        outputs = [];
%        rout = 100*ones(1,208);
% end

J = rout'*rout;
