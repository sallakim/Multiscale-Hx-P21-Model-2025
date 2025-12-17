function init = initialconditions(pars,data,eta_Vtot) 

%{ 
    This function approximates steady-state initial conditions. 
    Inputs: 
    pars        - vector of parameters 
    data        - input data structure with data and global parameters 
    Outputs: 
    init        - vector of initial conditions 
%} 


%% Fixed parameters 

fixpars = data.fixpars; 

% Minimal pressures (kPa) 
P_LV_m = fixpars(5); 
P_SA_m = fixpars(6); 
P_SV_m = fixpars(7);
P_RV_m = fixpars(8); 
P_PA_m = fixpars(9); 
P_PV_m = fixpars(10); 

% Reference sarcomere length (cm) 
Lsref = fixpars(13); % 2.0e-4

% Minimal elastances (kPa cm^(-3)) 
E_LV_m = fixpars(11);
E_RV_m = fixpars(12); 

% Crossbridge State Variables 
P1_0_LV = fixpars(28); % 0th moment state A1, LV
P1_1_LV = fixpars(29); % 1st moment state A1, LV
P1_2_LV = fixpars(30); % 2nd moment state A1, LV
P2_0_LV = fixpars(31); % 0th moment state A2, LV
P2_1_LV = fixpars(32); % 1st moment state A2, LV
P2_2_LV = fixpars(33); % 2nd moment state A2, LV
P3_0_LV = fixpars(34); % 0th moment state A3, LV
P3_1_LV = fixpars(35); % 1st moment state A3, LV
P3_2_LV = fixpars(36); % 2nd moment state A3, LV
N_LV = fixpars(37);
U_NR_LV = fixpars(38);
P1_0_SEP = fixpars(39); % 0th moment state A1, LV
P1_1_SEP = fixpars(40); % 1st moment state A1, LV
P1_2_SEP = fixpars(41); % 2nd moment state A1, LV
P2_0_SEP = fixpars(42); % 0th moment state A2, LV
P2_1_SEP= fixpars(43); % 1st moment state A2, LV
P2_2_SEP = fixpars(44); % 2nd moment state A2, LV
P3_0_SEP = fixpars(45); % 0th moment state A3, LV
P3_1_SEP = fixpars(46); % 1st moment state A3, LV
P3_2_SEP = fixpars(47); % 2nd moment state A3, LV
N_SEP = fixpars(48);
U_NR_SEP = fixpars(49);
P1_0_RV = fixpars(50); % 0th moment state A1, LV
P1_1_RV = fixpars(51); % 1st moment state A1, LV
P1_2_RV = fixpars(52); % 2nd moment state A1, LV
P2_0_RV = fixpars(53); % 0th moment state A2, LV
P2_1_RV = fixpars(54); % 1st moment state A2, LV
P2_2_RV = fixpars(55); % 2nd moment state A2, LV
P3_0_RV = fixpars(56); % 0th moment state A3, LV
P3_1_RV = fixpars(57); % 1st moment state A3, LV
P3_2_RV = fixpars(58); % 2nd moment state A3, LV
N_RV = fixpars(59);
U_NR_RV = fixpars(60);

%% Adjustable 
C_SA = pars(1); 
C_SV = pars(2); 
C_PA = pars(3); 
C_PV = pars(4); 

%% Initial conditions 

% Displacements (cm)                 1 - 4 
xm_LV_0  = data.deformation.xm_LV_0;
xm_SEP_0 = data.deformation.xm_SEP_0;
xm_RV_0  = data.deformation.xm_RV_0;
ym_0     = data.deformation.ym_0; 

% Sarcomere lengths (cm)             5 - 7 
Lsc_LV_0  = Lsref;
Lsc_SEP_0 = Lsref;
Lsc_RV_0  = Lsref;

% Volumes (cm^3)                     8 - 15 
V_LV_0 = eta_Vtot*P_LV_m / E_LV_m;
V_RV_0 = eta_Vtot*P_RV_m / E_RV_m; 
V_SV_0 = eta_Vtot*C_SV * P_SV_m;  
V_PV_0 = eta_Vtot*C_PV * P_PV_m;   
V_SA_0 = eta_Vtot*C_SA * P_SA_m;  
V_PA_0 = eta_Vtot*C_PA * P_PA_m;  

% Create initial conditions vector 
init = [xm_LV_0 ;xm_SEP_0 ;xm_RV_0 ;ym_0;           % 1-4
    Lsc_LV_0; Lsc_SEP_0; Lsc_RV_0;                  % 5-7
    V_LV_0; V_RV_0; V_SV_0; V_PV_0 ;V_SA_0 ;V_PA_0; % 8-13
    P1_0_LV; P1_1_LV; P1_2_LV ;P2_0_LV; P2_1_LV; P2_2_LV; P3_0_LV;P3_1_LV; P3_2_LV; N_LV; U_NR_LV;   % 14-24
    P1_0_SEP;P1_1_SEP;P1_2_SEP;P2_0_SEP;P2_1_SEP;P2_2_SEP;P3_0_SEP;P3_1_SEP;P3_2_SEP;N_SEP;U_NR_SEP; % 25-35
    P1_0_RV; P1_1_RV; P1_2_RV; P2_0_RV; P2_1_RV; P2_2_RV; P3_0_RV; P3_1_RV; P3_2_RV; N_RV; U_NR_RV]; % 36-46

% Solve system of Triseg equations to determine consistent initialization
% for DAE 
x0   = log(init(1:4)); 
opts = optimoptions('fsolve','Display','off',...
    'MaxFunctionEvaluations',2e3); 
xopt = fsolve(@(x) triseg(x,pars,data,init),x0,opts); 
init(1:4) = exp(xopt(1:4)); 

end 