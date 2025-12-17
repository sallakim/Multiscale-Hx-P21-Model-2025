function [adjpars,UB,LB,data] = parameters_Hx(data)

 %{ 
    Assignment and/or nominal calculation of all model parameters for Hx.. 
    
    Inputs: 
    data        - input data structure with data and global parameters 
    Outputs: 
    adjpars     - vector of adjustable parameters 
    UB          - vector of parameter upper bounds 
    LB          - vector of parameter lower bounds 
    data        - output data structure with new field assignments 
 %} 

%% Scaled for P21
eta_ka = 0.8;

%% Hx scaling 
eta_kstiff_RV    = 1.7; 
eta_k_off_RV     = 0.6;
eta_sigma_pas_RV = 2; 

%% Load in values from data structure 
% Blood pressures 
SPbar = data.SPbar; 
DPbar = data.DPbar; 

% Total blood volume (cm^3) 
Vtot  = data.Vtot; 

% Cardiac output (cm^3 s^(-1))
CO    = data.CO;

% End-diastolic and end-systolic pressures (kPa) and volumes (cm^3) 
ESP_LV = data.ESP_LV; 
EDP_LV = data.EDP_LV; 
EDV_LV = data.EDV_LV; 

ESP_RV = data.ESP_RV; 
EDP_RV = data.EDP_RV;   
EDV_RV = data.EDV_RV;   

% Wall volumes (mL) 
Wall_Volume_LV_and_SEP = data.Wall_Volume_LV_and_SEP;
Wall_Volume_RV         = data.Wall_Volume_RV;

%% Pericardium 

Vh0 = 1.25*(EDV_LV + EDV_RV); 
s = 10; 

%% Blood volume distribution 
% Fraction of unstressed volume in the compartments 
bvd_SA = .7;        bvd_PA = .4;
bvd_SV = .9;        bvd_PV = .9;

% bvd_SA = data.bvd_SA;   bvd_PA = data.bvd_PA; 
% bvd_SV = data.bvd_SV;   bvd_PV = data.bvd_PV; 

% Total blood volume distribution 
d_LV = .025;        d_RV = .025;
d_SA = .1;          d_PA = .05; 
d_SV = .45;         d_PV = .35;

d_sum_LV = d_PV + d_LV;
d_sum_RV = d_SV + d_RV;

d_LV = EDV_LV / Vtot; 
d_RV = EDV_RV / Vtot;
d_PV = d_sum_LV - d_LV;
d_SV = d_sum_RV - d_RV;

%% Calculate volumes (cm^3)

% Total chamber volumes: bld volume distribution fraction x total volume 
V_LV_0 = d_LV*Vtot;      V_RV_0 = d_RV*Vtot; 
V_SA_0 = d_SA*Vtot;      V_PA_0 = d_PA*Vtot; 
V_SV_0 = d_SV*Vtot;      V_PV_0 = d_PV*Vtot;

% Unstressed volumes
V_SA_u = V_SA_0*bvd_SA;   V_PA_u = V_PA_0*bvd_PA; 
V_SV_u = V_SV_0*bvd_SV;   V_PV_u = V_PV_0*bvd_PV; 

% Stressed volumes
V_SA_s = V_SA_0 - V_SA_u;  V_PA_s = V_PA_0 - V_PA_u;  
V_SV_s = V_SV_0 - V_SV_u;  V_PV_s = V_PV_0 - V_PV_u; 

%% Pressures (kPa)
% Max/min pressure ratios from Boron book and data (mmHg converted to kPa)
% LV                    % SA                                    % SV
P_LV_M  = ESP_LV;       P_SA_M   = SPbar;                       P_SV_M   = 8 / 7.5;
                        P_SA_bar = DPbar+1/3*(SPbar-DPbar);     P_SV_bar = 3 / 7.5; 
P_LV_m  = EDP_LV;       P_SA_m   = DPbar;                       P_SV_m   = 2 / 7.5;
                                            
% RV                    % PA                                    % PV
P_RV_M = ESP_RV;        P_PA_M   = ESP_RV/1.05;                 P_PV_M   = 8  / 7.5; 
                                                                P_PV_bar = 3 / 7.5; 
P_RV_m = EDP_RV;        P_PA_m   = P_PA_M - 15/7.5;             P_PV_m   = 2  / 7.5; 

%% Calculate elastances (kPa cm^(-3)) and compliances (cm^3 kPa^(-1))
 
% LV and RV minimal elastance parameters 
E_LV_m = P_LV_m / V_LV_0; 
E_RV_m = P_RV_m / V_RV_0; 

% Compliances
C_SA = V_SA_s/P_SA_M;
C_SV = 3;
C_PA = V_PA_s/P_PA_M; 
C_PV = V_PV_s/P_PV_M; 

%% Calculate resistances (kPa s cm^(-3))    
R_SA = (P_SA_M - P_SV_bar)/CO;
R_PA = (P_PA_M - P_PV_bar)/CO;

% Transmural resistances 
R_tSA =0; % 5.3432;%0.08 / 7.5e-0; % orig 0.05
R_tPA =0; % 1.3451;%0.02 / 7.5e-0; % orginal tPA = tSA

% Valves (vlv)
R_m = .1; 
R_a = 0.0333;
R_t = .1; 
R_p = 0.0333;

%% Heart model parameters 
% Conversion factor
mmHg_to_kPa = 1/7.5;

% Sarcomere geometry parameters (µm)
Lsref   = 2.0;  % Reference sarcomere length 
Lsc0    = 1.51; % Contractile element length 
L_thick = 1.67; % Length of thick filament
L_hbare = 0.10; % Length of bare region of thick filament
L_thin  = 1.20; % Length of thin filament

% Passive force parameters (mmHg/µm converted to kPa/µm)
k_passive_LV = 50 *mmHg_to_kPa;                   % for mean SHAM rat and TAC rat 1
k_passive_RV = eta_sigma_pas_RV * 50 *mmHg_to_kPa; 
gamma   = 7;                                      % steepness constant (unitless) 

% Active force parameters (kPa/µm)
kstiff1_LV = 1.4 * ( 1.7351e+03);
kstiff1_RV = eta_kstiff_RV * 1.4 * ( 1.7351e+03); 
kstiff2_LV = 1.4 * ( 4.5545e+04);     
kstiff2_RV = eta_kstiff_RV * 1.4 * ( 4.5545e+04); 
deltaR  = 0.010; % Crossbridge strain associated with ratcheting deformation (µm)

% Viscous force parameters (kPa s/µm)
eta    = 1*mmHg_to_kPa; % visconsity

% Series-elastic force parameters (kPa/µm)
Kse    = 50000*mmHg_to_kPa; % series element elastance

%% Calculate patient-specific reference midwall surface area (Amref) for LV, SEP, and RV (cm)

% Ventricular inner chamber radius 
r_LV_and_SEP = (EDV_LV * 3 / (4* pi))^(1/3); 
r_RV         = (EDV_RV * 3 / (4* pi))^(1/3); 

% Wall thickness (cm) 
if exist('Wall_Volume_RV','var') % determine wall thickness from inner chamber radius, EDV, and wall volume
    h_LV_and_SEP = ( (Wall_Volume_LV_and_SEP + EDV_LV) * 3/(4*pi)) ^ (1/3)-r_LV_and_SEP;
else
    h_LV_and_SEP = 0.3;
end

h_RV = h_LV_and_SEP/2; 

% Ventricle midwall radius (chamber radius (r) + 1/2 wall thickness (h)) where h is 
r_m_LV_and_SEP = r_LV_and_SEP + h_LV_and_SEP/2; % m 
r_m_RV         = r_RV + h_RV/2; % m 

% Outer radius 
r_o_LV_and_SEP = r_LV_and_SEP + h_LV_and_SEP; % m 
r_o_RV         = r_RV + h_RV; % m 

% Midwall reference surface area 
Amref_LV_and_SEP = 4 * pi * (r_m_LV_and_SEP)^2; 
Am_RV            = 4 * pi * (r_m_RV)^2;

Amref_LV  = Amref_LV_and_SEP * 2/3; % Assume LV is 2/3 of LV+SEP 
Amref_SEP = Amref_LV_and_SEP * 1/3; % Assume SEP is 1/3 of LV+SEP
Amref_RV  = Am_RV;

%% Calculate patient-specific midwall volume (Vw) for LV, SEP, and RV 

% Outer Ventricle volume 
Vw_chamber_LV_and_SEP = 4/3 * pi * r_o_LV_and_SEP^3;  
Vw_chamber_RV         = 4/3 * pi * r_o_RV^3; 

% Ventricular wall volume = weight (g) / 1.055(g/mL) 
if exist('Wall_Volume_RV','var')
    Vw_RV = Wall_Volume_RV; 
else
    Vw_RV = Vw_chamber_RV - EDV_RV;
end

if exist('Wall_Volume_LV_and_SEP', 'var')
    Vw_LV  = 2/3 * Wall_Volume_LV_and_SEP;
    Vw_SEP = 1/3 * Wall_Volume_LV_and_SEP;
else
    Vw_LV_and_SEP = Vw_chamber_LV_and_SEP - EDV_LV;
    Vw_LV =  Vw_LV_and_SEP * 2/3;
    Vw_SEP = Vw_LV_and_SEP * 1/3;
end

%% Approximations for initial displacements and Amref_rv in end-diastole 

% Initialize diastolic displacement values 
xm_LV_d_0  = 0.32; 
xm_SEP_d_0 = 0.13; 
xm_RV_d_0  = 0.38; 
ym_d_0     = 0.19; 

x0 = [xm_LV_d_0; 
    xm_SEP_d_0; 
    xm_RV_d_0;
    ym_d_0; 
    Amref_RV; 
    ]; 
x0 = log(x0); % log-scale the initial values 

% Inputs for calculating displacements 
Vw    = [Vw_LV,Vw_SEP,Vw_RV]; 
Amref = [Amref_LV,Amref_SEP]; 

% Assume end-diastolic sarcomere length (µm) 
SL_d    = 2.0;  

opts = optimoptions('fsolve','Display','none',...
    'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
[fnew0,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,SL_d,EDV_LV,0),x0,opts); 
fnew0 = exp(fnew0);

% Outputs / Diastolic displacements (cm)  
xm_LV_d  = fnew0(1);
xm_SEP_d = fnew0(2);
xm_RV_d  = fnew0(3);
ym_d     = fnew0(4);
Amref_RV = fnew0(5); 

% Set initial conditions for displacements to be used in the
% initialconditions.m script 
deformation.xm_LV_0  = xm_LV_d; 
deformation.xm_SEP_0 = xm_SEP_d; 
deformation.xm_RV_0  = xm_RV_d; 
deformation.ym_0     = ym_d; 

data.deformation = deformation; 

%% Calcium Moments 

P1_0_LV = 0; % 0th moment state A1, LV
P1_1_LV = 0; % 1st moment state A1, LV
P1_2_LV = 0; % 2nd moment state A1, LV
P2_0_LV = 0; % 0th moment state A2, LV
P2_1_LV = 0; % 1st moment state A2, LV
P2_2_LV = 0; % 2nd moment state A2, LV
P3_0_LV = 0; % 0th moment state A3, LV
P3_1_LV = 0; % 1st moment state A3, LV
P3_2_LV = 0; % 2nd moment state A3, LV
N_LV = 1;
U_NR_LV = 0;
P1_0_SEP = 0; % 0th moment state A1, LV
P1_1_SEP = 0; % 1st moment state A1, LV
P1_2_SEP = 0; % 2nd moment state A1, LV
P2_0_SEP = 0; % 0th moment state A2, LV
P2_1_SEP = 0; % 1st moment state A2, LV
P2_2_SEP = 0; % 2nd moment state A2, LV
P3_0_SEP = 0; % 0th moment state A3, LV
P3_1_SEP = 0; % 1st moment state A3, LV
P3_2_SEP = 0; % 2nd moment state A3, LV
N_SEP = 1;
U_NR_SEP = 0;
P1_0_RV = 0; % 0th moment state A1, LV
P1_1_RV = 0; % 1st moment state A1, LV
P1_2_RV = 0; % 2nd moment state A1, LV
P2_0_RV = 0; % 0th moment state A2, LV
P2_1_RV = 0; % 1st moment state A2, LV
P2_2_RV = 0; % 2nd moment state A2, LV
P3_0_RV = 0; % 0th moment state A3, LV
P3_1_RV = 0; % 1st moment state A3, LV
P3_2_RV = 0; % 2nd moment state A3, LV
N_RV = 1;
U_NR_RV = 0;

%% Crossbridge 

% Permissible (P) to loosely-bound (A1) 
ka      = eta_ka*559.5568; % myosin-actin attach rate constant, 1/sec
kd      = 304.6708; % myosin-actin detach rate constant, 1/sec

% Loosely-bound (A1, A1') to strongly-bound (A3) 
K_Pi    = 4.00;     % Pi unattachment rate 
k1      = 112.3727; % transition A1 to A2 rate constant, 1/sec
km1     = 21.296;   % transition A2 to A1 rate constant, 1/sec
alpha1  = 10.0;     % Stretch sensing parameter for k1 and km1, 1/µm 

% Strongly-bound (A2) to post-ratcheted (A3) 
k2      = 811.72;   % transition A2 to A3 rate constant, 1/sec
km2     = 43.25;    % transition A3 to A2 rate constant, 1/sec
alpha2  = 9.1;      % Stretch sensing parameter for k2 and km2, 1/µm  

% Post-ratcheted states (A3, A3', A3'') to permissible (P)
K_D     = 0.194;    % ADP unattachment rate, µm
K_T     = 0.4897;   % ATP unattachment rate, µm
k3      = 144.5586; % transition A3 to P rate constant, 1/sec
alpha3  = 0.1*59.3; % Stretch sensing parameter for k3, 1/µm 
s3      = 9.9e-3;   % Strain at which k3 is minimum, µm

% Calcium activation, permissible (P) to non-permissible (N) 
k_on     = 101.1850;   % Campbell et al, Biophysical Journal 2018
k_off_LV = 723.8520;   % manually tuned parameter!
k_off_RV = eta_k_off_RV*723.8520;     % manually tuned parameter!
K_coop   = 9.6846;     % Campbell et al, Biophysical Journal 2018,

% Super relaxed state to non relaxed state
ksr    = 15;      % s^(-2)
kmsr   = 50.032;  % s^(-2)
kforce = 1.169;   % N^(-1) m^(-2)

% Calcium transient timings 
k_TS = 0.1; % Rise
k_TR = 0.3; % Fall 

%% Outputs

parameters = [
    ka; kd;                             % 1-2                                          
    K_Pi; k1; km1; alpha1;              % 3-6           
    k2; km2; alpha2;                    % 7-9          
    K_D; K_T; k3; alpha3; s3;           % 10-14                
    k_on; k_off_LV; k_off_RV; K_coop;   % 15-18                  
    ksr; kmsr; kforce;                  % 19-21
    k_TS; k_TR;                         % 22-23

    Lsref; Lsc0; L_thick; L_hbare; L_thin;           % 24-28   
    k_passive_LV; k_passive_RV; gamma;               % 29-31
    kstiff1_LV; kstiff1_RV; kstiff2_LV; kstiff2_RV;  % 32-35
    deltaR; eta; Kse;                                % 36-38      

    Vw_LV; Vw_SEP; Vw_RV;           % 39-41
    Amref_LV; Amref_SEP; Amref_RV;  % 42-44

    V_SA_u; V_SV_u; V_PA_u; V_PV_u;                 % 45-48                 
    P_LV_m; P_SA_m; P_SV_m; P_RV_m; P_PA_m; P_PV_m; % 49-54
    E_LV_m; E_RV_m;                                 % 55-56
    C_SA; C_SV; C_PA; C_PV;                         % 57-60
    R_SA; R_PA;                                     % 61-62
    R_m; R_a; R_t; R_p;                             % 63-66

    P1_0_LV; P1_1_LV; P1_2_LV; P2_0_LV; P2_1_LV; P2_2_LV; P3_0_LV; P3_1_LV; P3_2_LV; N_LV; U_NR_LV;             % 67-77      
    P1_0_SEP;P1_1_SEP; P1_2_SEP; P2_0_SEP; P2_1_SEP; P2_2_SEP; P3_0_SEP; P3_1_SEP; P3_2_SEP; N_SEP; U_NR_SEP;   % 78-88
    P1_0_RV; P1_1_RV; P1_2_RV; P2_0_RV; P2_1_RV; P2_2_RV; P3_0_RV; P3_1_RV; P3_2_RV; N_RV; U_NR_RV;             % 89-99
    ];

adjpars = [parameters(57:66);parameters(42:44);parameters(39:41);parameters(22:23)];
  
% adjpars = [C_SA; C_SV; C_PA; C_PV;      % 1-4
%     R_SA; R_PA;                         % 5-6
%     R_m; R_a; R_t; R_p;                 % 7-10
% 
%     Amref_LV; Amref_SEP; Amref_RV;      % 11-13     
%     Vw_LV; Vw_SEP; Vw_RV;               % 14-16
% 
%     k_TS; k_TR;                         % 17-18
%     ]; 

UB = adjpars*10; 
LB = adjpars/10;

UB(11:16) = adjpars(11:16)*1.2;
LB(11:16) = adjpars(11:16)*0.8;

UB(17) = .001;
LB(17) = 1;

UB(18) = .001;
LB(18) = 1;

adjpars = log(adjpars); 
UB = log(UB);
LB = log(LB);

fixpars = [parameters(45:56);parameters(24:38);parameters(67:99);parameters(1:21);Vh0;s;R_tSA;R_tPA];

% fixpars = [V_SA_u; V_SV_u; V_PA_u; V_PV_u;              % 1-4       
%     P_LV_m; P_SA_m; P_SV_m; P_RV_m; P_PA_m; P_PV_m;     % 5-10
%     E_LV_m; E_RV_m;                                     % 11-12
%     Lsref; Lsc0; L_thick; L_hbare; L_thin;              % 13-17
%     k_passive_LV; k_passive_RV; gamma;                  % 18-20
%     kstiff1_LV; kstiff1_RV; kstiff2_LV; kstiff2_RV; deltaR; % 21-25
%     eta; Kse;                                           % 26-27
%     P1_0_LV; P1_1_LV; P1_2_LV; P2_0_LV; P2_1_LV; P2_2_LV; P3_0_LV; P3_1_LV; P3_2_LV; N_LV; U_NR_LV;           % 28-38
%     P1_0_SEP;P1_1_SEP; P1_2_SEP; P2_0_SEP; P2_1_SEP; P2_2_SEP; P3_0_SEP; P3_1_SEP; P3_2_SEP; N_SEP; U_NR_SEP; % 39-49
%     P1_0_RV; P1_1_RV; P1_2_RV; P2_0_RV; P2_1_RV; P2_2_RV; P3_0_RV; P3_1_RV; P3_2_RV; N_RV; U_NR_RV;           % 50-60 
%     ka; kd;                                             % 61-62
%     K_Pi; k1; km1; alpha1;                              % 63-66
%     k2; km2; alpha2;                                    % 67-69
%     K_D; K_T; k3; alpha3; s3;                           % 70-74
%     k_on; k_off_LV; k_off_RV; K_coop;                   % 75-78
%     ksr; kmsr; kforce;                                  % 79-81
%     Vh0; s;                                             % 82-83
%     R_tSA; R_tPA];                                      % 84-85

data.fixpars = fixpars; 

end 