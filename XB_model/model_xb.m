function [dxdt,outputs] = model_xb(t,init,pars,data)

if isfield(data,'eta_sigma_pas') == 1 
    eta_sigma_pas = 2;
else
    eta_sigma_pas = 1; 
end


%% Data 
MgATP_LV = data.MgATP_cyto;
MgATP_SEP = data.MgATP_cyto;
MgATP_RV = data.MgATP_cyto;
MgADP_LV = data.MgADP_cyto;
MgADP_SEP = data.MgADP_cyto;
MgADP_RV = data.MgADP_cyto;
Pi_LV = data.Pi_cyto;
Pi_SEP = data.Pi_cyto;
Pi_RV = data.Pi_cyto;
T = data.T;

% Calcium Concentration 
if isfield(data,'Ca_i') == 1 
    Ca_i = data.Ca_i; 
else
    Ca_amplitude = 2;
    Ca_diastole = 0.1610;
    tc = mod(t,T);
    TS_v = 0.1 * T; 
    TR_v = 0.3 * T; 
    if tc >= 0 && tc < TS_v
        Ca_i = Ca_diastole + Ca_amplitude*0.5*(1 - cos(pi*(tc)/TS_v)); 
    elseif tc >= TS_v && tc < TR_v + TS_v 
        Ca_i = Ca_diastole + Ca_amplitude*0.5*(1 + cos(pi*(tc - TS_v)/TR_v)); 
    else
        Ca_i = Ca_diastole; 
    end 
end


% Sarcomere length (um)
Ls_LV  = data.Ls; 
Ls_SEP = data.Ls; 
Ls_RV  = data.Ls; 

%% init 
Lsc_LV = init(1); 
Lsc_SEP = init(2); 
Lsc_RV = init(3);

% Calcium moments 
P1_0_LV = init(4); % 0th moment state A1, LV
P1_1_LV = init(5); % 1st moment state A1, LV
P1_2_LV = init(6); % 2nd moment state A1, LV
P2_0_LV = init(7); % 0th moment state A2, LV
P2_1_LV = init(8); % 1st moment state A2, LV
P2_2_LV = init(9); % 2nd moment state A2, LV
P3_0_LV = init(10); % 0th moment state A3, LV
P3_1_LV = init(11); % 1st moment state A3, LV
P3_2_LV = init(12); % 2nd moment state A3, LV
N_LV = init(13);
U_NR_LV = init(14);
P1_0_SEP = init(15); % 0th moment state A1, LV
P1_1_SEP = init(16); % 1st moment state A1, LV
P1_2_SEP = init(17); % 2nd moment state A1, LV
P2_0_SEP = init(18); % 0th moment state A2, LV
P2_1_SEP = init(19);  % 1st moment state A2, LV
P2_2_SEP = init(20); % 2nd moment state A2, LV
P3_0_SEP = init(21); % 0th moment state A3, LV
P3_1_SEP = init(22); % 1st moment state A3, LV
P3_2_SEP = init(23); % 2nd moment state A3, LV
N_SEP = init(24);
U_NR_SEP = init(25);
P1_0_RV = init(26); % 0th moment state A1, LV
P1_1_RV = init(27); % 1st moment state A1, LV
P1_2_RV = init(28); % 2nd moment state A1, LV
P2_0_RV = init(29); % 0th moment state A2, LV
P2_1_RV = init(30); % 1st moment state A2, LV
P2_2_RV = init(31); % 2nd moment state A2, LV
P3_0_RV = init(32); % 0th moment state A3, LV
P3_1_RV = init(33); % 1st moment state A3, LV
P3_2_RV = init(34); % 2nd moment state A3, LV
N_RV = init(35);
U_NR_RV = init(36);

%% Parameters 
% Collagen Froce 
SLcollagen   = pars(1); % threshold for collagen activation,(um)
PConcollagen = pars(2); % contriubtion of collagen (unitless)
PExpcollagen = pars(3); % expresion of collagen (unitless)

% Passive Stress
eta        = pars(4); % visconsity, mmHg s /micron
L_rest_pas = pars(5); % (um)  

kstiff1 = pars(6);   % kPa/um (9/5 BM)
kstiff2 = pars(7);   % kPa/um (9/5 BM)
k_passive = pars(8); % kPa / um % for mean SHAM rat and TAC rat 1

% Sarcomere Geometry Parameters (um)
L_thick = pars(9);  % Length of thick filament
L_hbare = pars(10); % Length of bare region of thick filament
L_thin  = pars(11); % Length of thin filament
deltaR  = pars(12); 

alpha1  = pars(13); % Stretch sensing parameter for k1 and k?1, 1/um 
alpha2  = pars(14); % Stretch sensing parameter for k2 and k?2, 1/um  
alpha3  = pars(15); % Stretch sensing parameter for k3, 1/um 
s3      = pars(16); % Strain at which k3 is minimum, um
km1     = pars(17); % transition A2 to A1 rate constant, 1/sec
km2     = pars(18); % transition A3 to A2 rate constant, 1/sec
k1      = pars(19); % transition A1 to A2 rate constant, 1/sec
k2      = pars(20); % transition A2 to A3 rate constant, 1/sec
k3      = pars(21); % transition A3 to P rate constant, 1/sec
ka      = pars(22); % myosin-actin attach rate constant, 1/sec
kd      = pars(23); % myosin-actin detach rate constant, 1/sec
ksr     = pars(24);
kforce  = pars(25); % dived by kPa to mmHg conversion rate
kmsr    = pars(26); % fit to the invitro data
k_on    = pars(27); % Campbell et al, Biophysical Journal 2018
K_coop  = pars(28); % Campbell et al, Biophysical Journal 2018,
Kse     = pars(29); % series element elastance, mmHg/micron to kPa/micron 
K_Pi    = pars(30); 
K_T     = pars(31); 
K_D     = pars(32); % Used the values from Tewari etal JMCC (9/5 BM)
k_off   = pars(33);

%% Correcting rate constants 
kd_LV  = kd*(Pi_LV/K_Pi)/(1.0 + Pi_LV/K_Pi);
k1_LV  = k1/(1.0 + Pi_LV/K_Pi);
km2_LV = km2*(MgADP_LV/K_D)/(1.0 + MgADP_LV/K_D + MgATP_LV/K_T);
k3_LV  = k3*(MgATP_LV/K_T)/(1.0 + MgATP_LV/K_T + MgADP_LV/K_D);

kd_SEP  = kd*(Pi_SEP/K_Pi)/(1.0 + Pi_SEP/K_Pi);
k1_SEP  = k1/(1.0 + Pi_SEP/K_Pi);
km2_SEP = km2*(MgADP_SEP/K_D)/(1.0 + MgADP_SEP/K_D + MgATP_SEP/K_T);
k3_SEP  = k3*(MgATP_SEP/K_T)/(1.0 + MgATP_SEP/K_T + MgADP_SEP/K_D);

kd_RV  = kd*(Pi_RV/K_Pi)/(1.0 + Pi_RV/K_Pi);
k1_RV  = k1/(1.0 + Pi_RV/K_Pi);
km2_RV = km2*(MgADP_RV/K_D)/(1.0 + MgADP_RV/K_D + MgATP_RV/K_T);
k3_RV  = k3*(MgATP_RV/K_T)/(1.0 + MgATP_RV/K_T + MgADP_RV/K_D);

%% XB Stress 
% Sarcomere geometry (um) 
sovr_ze = min(L_thick*0.5, Lsc_LV*0.5);
sovr_cle = max(Lsc_LV*0.5 - (Lsc_LV-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_LV = L_sovr*2/(L_thick - L_hbare);

sovr_ze = min(L_thick*0.5, Lsc_SEP*0.5);
sovr_cle = max(Lsc_SEP*0.5 - (Lsc_SEP-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_SEP = L_sovr*2/(L_thick - L_hbare);

sovr_ze = min(L_thick*0.5, Lsc_RV*0.5);
sovr_cle = max(Lsc_RV*0.5 - (Lsc_RV-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_RV = L_sovr*2/(L_thick - L_hbare); % unitless 

% Collagen force
% sigma_collagen_LV  = PConcollagen*(exp(PExpcollagen*(Ls_LV - SLcollagen)) - 1).*(Ls_LV > SLcollagen);
% sigma_collagen_SEP = PConcollagen*(exp(PExpcollagen*(Ls_SEP - SLcollagen)) - 1).*(Ls_SEP > SLcollagen);
% sigma_collagen_RV  = PConcollagen*(exp(PExpcollagen*(Ls_RV - SLcollagen)) - 1).*(Ls_RV > SLcollagen);

% Passive Stress (kPa) 
% sigma_pas_LV  = k_passive*(Ls_LV/2-L_rest_pas)  ;
% sigma_pas_SEP = k_passive*(Ls_SEP/2-L_rest_pas) ;
% sigma_pas_RV  = eta_sigma_pas*k_passive*(Ls_RV/2-L_rest_pas)  ;
% 
% sigma_pas_LV  = sigma_pas_LV  + sigma_collagen_LV ;
% sigma_pas_SEP = sigma_pas_SEP + sigma_collagen_SEP;
% sigma_pas_RV  = sigma_pas_RV  + sigma_collagen_RV;

gamma = 7;
Lsc0    = 1.51;
sigma_pas_LV  =  k_passive * ((Ls_LV - Lsc0))^gamma; 
sigma_pas_SEP =  k_passive * ((Ls_SEP - Lsc0))^gamma; 
sigma_pas_RV  =  eta_sigma_pas * k_passive * ((Ls_RV - Lsc0))^gamma; 

% Active stress (kPa) from XB 
sigma_act_LV  = N_overlap_LV*(kstiff2*deltaR*(P3_0_LV) + kstiff1*(P2_1_LV + P3_1_LV)); % mmHg * normalised force
sigma_act_SEP = N_overlap_SEP*(kstiff2*deltaR*(P3_0_SEP) + kstiff1*(P2_1_SEP+P3_1_SEP)); % mmHg * normalised force
sigma_act_RV  = N_overlap_RV*(kstiff2*deltaR*(P3_0_RV) + kstiff1*(P2_1_RV+P3_1_RV)); % mmHg * normalised force

% Total stress (kPa)
sigma_LV  = -Kse*(Lsc_LV - Ls_LV);
sigma_SEP = -Kse*(Lsc_SEP - Ls_SEP);
sigma_RV  = -Kse*(Lsc_RV - Ls_RV);

%% Myofiber Mechanics (Lsc_LV) 
% Calculations for stretch-senstive rates    
f_alpha1o_LV = (P1_0_LV - alpha1*P1_1_LV + 0.5*(alpha1*alpha1)*P1_2_LV);
f_alpha1i_LV = (P1_1_LV - alpha1*P1_2_LV);

f_alpha0o_LV = (P2_0_LV + alpha1*P2_1_LV + 0.5*alpha1*alpha1*P2_2_LV);
f_alpha0i_LV = (P2_1_LV + alpha1*P2_2_LV);

f_alpha2o_LV = (P2_0_LV - alpha2*P2_1_LV + 0.5*(alpha2*alpha2)*P2_2_LV);
f_alpha2i_LV = (P2_1_LV - alpha2*P2_2_LV);

f_alpha3o_LV = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV));
f_alpha3i_LV = (P3_1_LV + alpha3*(s3*s3*P3_1_LV + 2.0*s3*P3_2_LV));

%% Myofiber Mechanics (Lsc_SEP) 
% Calculations for stretch-senstive rates    
f_alpha1o_SEP = (P1_0_SEP - alpha1*P1_1_SEP + 0.5*(alpha1*alpha1)*P1_2_SEP);
f_alpha1i_SEP = (P1_1_SEP - alpha1*P1_2_SEP);

f_alpha0o_SEP = (P2_0_SEP + alpha1*P2_1_SEP + 0.5*alpha1*alpha1*P2_2_SEP);
f_alpha0i_SEP = (P2_1_SEP + alpha1*P2_2_SEP);

f_alpha2o_SEP = (P2_0_SEP - alpha2*P2_1_SEP + 0.5*(alpha2*alpha2)*P2_2_SEP);
f_alpha2i_SEP = (P2_1_SEP - alpha2*P2_2_SEP);

f_alpha3o_SEP = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP));
f_alpha3i_SEP = (P3_1_SEP + alpha3*(s3*s3*P3_1_SEP + 2.0*s3*P3_2_SEP));

%% Myofiber Mechanics (Lsc_RV) 
% Calculations for stretch-senstive rates    
f_alpha1o_RV = (P1_0_RV - alpha1*P1_1_RV + 0.5*(alpha1*alpha1)*P1_2_RV);
f_alpha1i_RV = (P1_1_RV - alpha1*P1_2_RV);

f_alpha0o_RV = (P2_0_RV + alpha1*P2_1_RV + 0.5*alpha1*alpha1*P2_2_RV);
f_alpha0i_RV = (P2_1_RV + alpha1*P2_2_RV);

f_alpha2o_RV = (P2_0_RV - alpha2*P2_1_RV + 0.5*(alpha2*alpha2)*P2_2_RV);
f_alpha2i_RV = (P2_1_RV - alpha2*P2_2_RV);

f_alpha3o_RV = (P3_0_RV + alpha3*(s3*s3*P3_0_RV + 2.0*s3*P3_1_RV + P3_2_RV));
f_alpha3i_RV = (P3_1_RV + alpha3*(s3*s3*P3_1_RV + 2.0*s3*P3_2_RV));

%% XB 
% LV
P0_LV = 1.0 - N_LV - P1_0_LV - P2_0_LV - P3_0_LV; % DAB 10/8/2019

U_SR_LV = 1 - U_NR_LV;
Jon_LV = k_on*Ca_i*N_LV*(1 + K_coop*(1 - N_LV));
Joff_LV = k_off*P0_LV*(1 + K_coop*N_LV);
% SEP 
P0_SEP = 1 - N_SEP - P1_0_SEP - P2_0_SEP - P3_0_SEP;  % DAB 10/8/2019

U_SR_SEP = 1 - U_NR_SEP;
Jon_SEP = k_on*Ca_i*N_SEP*(1 + K_coop*(1 - N_SEP));
Joff_SEP = k_off*P0_SEP*(1 + K_coop*N_SEP);
% RV
P0_RV = 1.0 - N_RV - P1_0_RV - P2_0_RV - P3_0_RV;  % DAB 10/8/2019

U_SR_RV = 1.0 - U_NR_RV;
Jon_RV = k_on*Ca_i*N_RV*(1 + K_coop*(1 - N_RV));
Joff_RV = k_off*P0_RV*(1 + K_coop*N_RV);

%% dxdt
% (Myofiber Mechanics)
dLsc_LV = (sigma_LV - sigma_pas_LV - sigma_act_LV)/eta;
dLsc_SEP = (sigma_SEP - sigma_pas_SEP - sigma_act_SEP)/eta;
dLsc_RV = (sigma_RV - sigma_pas_RV - sigma_act_RV)/eta;

% (xb, LV) 
dP1_0_LV = ka*P0_LV*U_NR_LV*N_overlap_LV - kd_LV*P1_0_LV - k1_LV*f_alpha1o_LV + km1*f_alpha0o_LV; % DAB 10/8/2019
dP1_1_LV = dLsc_LV*P1_0_LV - kd_LV*P1_1_LV - k1_LV*f_alpha1i_LV + km1*f_alpha0i_LV;
dP1_2_LV = 2*dLsc_LV*P1_1_LV - kd_LV*P1_2_LV - k1_LV*P1_2_LV + km1*P2_2_LV;

dP2_0_LV = -km1*f_alpha0o_LV - k2*f_alpha2o_LV + km2_LV*P3_0_LV + k1_LV*f_alpha1o_LV;
dP2_1_LV = dLsc_LV*P2_0_LV - km1*f_alpha0i_LV - k2*f_alpha2i_LV + km2_LV*P3_1_LV + k1_LV*f_alpha1i_LV;
dP2_2_LV = 2*dLsc_LV*P2_1_LV - km1*P2_2_LV - k2*P2_2_LV + km2_LV*P3_2_LV + k1_LV*P1_2_LV;

dP3_0_LV = +k2*f_alpha2o_LV - km2_LV*P3_0_LV - k3_LV*f_alpha3o_LV;
dP3_1_LV = dLsc_LV*P3_0_LV + k2*f_alpha2i_LV - km2_LV*P3_1_LV - k3_LV*f_alpha3i_LV;
dP3_2_LV = 2*dLsc_LV*P3_1_LV + k2*P2_2_LV - km2_LV*P3_2_LV - k3_LV*P3_2_LV;

dN_LV    = - Jon_LV + Joff_LV; % dN_LV / dt
dU_NR_LV = ksr * (1 + kforce * sigma_act_LV) * U_SR_LV - kmsr*U_NR_LV ; % DAB 10/8

% (xb, SEP) 
dP1_0_SEP = ka*P0_SEP*U_NR_SEP*N_overlap_SEP - kd_SEP*P1_0_SEP - k1_SEP*f_alpha1o_SEP + km1*f_alpha0o_SEP; % DAB 10/8/2019
dP1_1_SEP = dLsc_SEP*P1_0_SEP - kd_SEP*P1_1_SEP - k1_SEP*f_alpha1i_SEP + km1*f_alpha0i_SEP;
dP1_2_SEP = 2*dLsc_SEP*P1_1_SEP - kd_SEP*P1_2_SEP - k1_SEP*P1_2_SEP + km1*P2_2_SEP;

dP2_0_SEP = -km1*f_alpha0o_SEP - k2*f_alpha2o_SEP + km2_SEP*P3_0_SEP + k1_SEP*f_alpha1o_SEP;
dP2_1_SEP = dLsc_SEP*P2_0_SEP - km1*f_alpha0i_SEP - k2*f_alpha2i_SEP + km2_SEP*P3_1_SEP + k1_SEP*f_alpha1i_SEP;
dP2_2_SEP = 2*dLsc_SEP*P2_1_SEP - km1*P2_2_SEP       - k2*P2_2_SEP + km2_SEP*P3_2_SEP + k1_SEP*P1_2_SEP;

dP3_0_SEP = +k2*f_alpha2o_SEP - km2_SEP*P3_0_SEP - k3_SEP*f_alpha3o_SEP;
dP3_1_SEP = dLsc_SEP*P3_0_SEP + k2*f_alpha2i_SEP - km2_SEP*P3_1_SEP - k3_SEP*f_alpha3i_SEP;
dP3_2_SEP = 2*dLsc_SEP*P3_1_SEP + k2*P2_2_SEP       - km2_SEP*P3_2_SEP - k3_SEP*P3_2_SEP;

dN_SEP    = - Jon_SEP + Joff_SEP; 
dU_NR_SEP = ksr*(1 + kforce * sigma_act_SEP) * U_SR_SEP - kmsr*U_NR_SEP ; % DAB 10/8

% (xb, RV) 
dP1_0_RV = ka*P0_RV*U_NR_RV*N_overlap_RV - kd_RV*P1_0_RV - k1_RV*f_alpha1o_RV + km1*f_alpha0o_RV;  % DAB 10/8/2019
dP1_1_RV = dLsc_RV*P1_0_RV - kd_RV*P1_1_RV - k1_RV*f_alpha1i_RV + km1*f_alpha0i_RV;
dP1_2_RV = 2*dLsc_RV*P1_1_RV - kd_RV*P1_2_RV - k1_RV*P1_2_RV + km1*P2_2_RV;

dP2_0_RV = -km1*f_alpha0o_RV - k2*f_alpha2o_RV + km2_RV*P3_0_RV + k1_RV*f_alpha1o_RV;
dP2_1_RV = dLsc_RV*P2_0_RV - km1*f_alpha0i_RV - k2*f_alpha2i_RV + km2_RV*P3_1_RV + k1_RV*f_alpha1i_RV;
dP2_2_RV = 2*dLsc_RV*P2_1_RV - km1*P2_2_RV       - k2*P2_2_RV + km2_RV*P3_2_RV + k1_RV*P1_2_RV;

dP3_0_RV = +k2*f_alpha2o_RV - km2_RV*P3_0_RV - k3_RV*f_alpha3o_RV;
dP3_1_RV = dLsc_RV*P3_0_RV + k2*f_alpha2i_RV - km2_RV*P3_1_RV - k3_RV*f_alpha3i_RV;
dP3_2_RV = 2*dLsc_RV*P3_1_RV + k2*P2_2_RV       - km2_RV*P3_2_RV - k3_RV*P3_2_RV;

dN_RV    = - Jon_RV + Joff_RV; 
dU_NR_RV = ksr*(1 + kforce * sigma_act_RV) * U_SR_RV - kmsr * U_NR_RV ; % DAB 10/8

dxdt = [
    dLsc_LV; dLsc_SEP; dLsc_RV;                 %   
    dP1_0_LV; dP1_1_LV; dP1_2_LV;               % 
    dP2_0_LV; dP2_1_LV; dP2_2_LV;               % 
    dP3_0_LV; dP3_1_LV; dP3_2_LV;               % 
    dN_LV; dU_NR_LV;                            % 
    dP1_0_SEP; dP1_1_SEP; dP1_2_SEP;            % 
    dP2_0_SEP; dP2_1_SEP; dP2_2_SEP;            % 
    dP3_0_SEP; dP3_1_SEP; dP3_2_SEP;            % 
    dN_SEP; dU_NR_SEP;                          % 
    dP1_0_RV; dP1_1_RV; dP1_2_RV;               % 
    dP2_0_RV; dP2_1_RV; dP2_2_RV;               % 
    dP3_0_RV; dP3_1_RV; dP3_2_RV;               % 
    dN_RV; dU_NR_RV;                            % 
    ];  

%% Outputs 
outputs = [sigma_act_LV, sigma_act_SEP, sigma_act_RV, sigma_LV, sigma_SEP, sigma_RV, sigma_pas_LV, sigma_pas_SEP, sigma_pas_RV];
