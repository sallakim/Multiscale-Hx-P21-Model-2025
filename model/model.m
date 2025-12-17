function [dxdt, outputs] = model(t,x,pars,data) 

%{ 
    This function contains the right-hand side of the in vivo model. 
    Inputs: 
    t       - time 
    x       - states 
    pars    - vector of adjustable parameter values 
    data    - input data structure with data and global parameters 
    Outputs: 
    dxdt    - vector of solved ODEs
    outputs - vector of solved auxiliary equations 
%}

%% Get Data 
fixpars = data.fixpars;

if isfield(data,'hx_flag') == 1 
    eta_ATP = 0.95;
    eta_ADP_Pi = 1.05;
else
    eta_ATP = 1;
    eta_ADP_Pi = 1;
end

HR = data.HR;
T  = data.T; 

MgATP_LV  = data.MgATP_cytoplasm;
MgATP_SEP = data.MgATP_cytoplasm;
MgATP_RV  = eta_ATP*data.MgATP_cytoplasm;

MgADP_LV  = data.MgADP_cytoplasm;
MgADP_SEP = data.MgADP_cytoplasm;
MgADP_RV  = eta_ADP_Pi*data.MgADP_cytoplasm;

Pi_LV  = data.Pi_cyto;
Pi_SEP = data.Pi_cyto;
Pi_RV  = eta_ADP_Pi*data.Pi_cyto;

%% Adjustable parameters

% Compliance (cm^3 kPa^(-1))
C_SA = pars(1); 
C_SV = pars(2); 
C_PA = pars(3); 
C_PV = pars(4); 

% Resistance (kPa s cm^(-3))
R_SA = pars(5); 
R_PA = pars(6); 

R_m_vlv = pars(7); 
R_a_vlv = pars(8); 
R_t_vlv = pars(9); 
R_p_vlv = pars(10); 

% Midwall reference surface area (cm^2)
Amref_LV  = pars(11); 
Amref_SEP = pars(12);
Amref_RV  = pars(13); 

% Midwall volume (cm^3)
Vw_LV  = pars(14); 
Vw_SEP = pars(15); 
Vw_RV  = pars(16); 

% Percentage of cardiac cycle 
k_TS = pars(17); % Beginning of cardiac cycle to maximal systole  
k_TR = pars(18); 

%% Fixed parameters 

% Sarcomere Geometry Parameters (um)
Lsref   = fixpars(13);
Lsc0    = fixpars(14); 
L_thick = fixpars(15); % Length of thick filament
L_hbare = fixpars(16); % Length of bare region of thick filament
L_thin  = fixpars(17); % Length of thin filament

% Passive 
k_passive_LV = fixpars(18); % kPa / um % for mean SHAM rat and TAC rat 1
k_passive_RV = fixpars(19); % kPa / um % for mean SHAM rat and TAC rat 1
gamma   = fixpars(20); 

% Active 
kstiff1_LV = fixpars(21);   % kPa/um (9/5 BM)
kstiff1_RV = fixpars(22);   % kPa/um (9/5 BM)
kstiff2_LV = fixpars(23);   % kPa/um (9/5 BM)
kstiff2_RV = fixpars(24);   % kPa/um (9/5 BM)
deltaR  = fixpars(25); 
eta       = fixpars(26); % visconsity, mmHg s /micron
Kse     = fixpars(27); % series element elastance, mmHg/micron to kPa/micron 

% Cycling rates 
ka      = fixpars(61); % myosin-actin attach rate constant, 1/sec
kd      = fixpars(62); % myosin-actin detach rate constant, 1/sec
K_Pi    = fixpars(63); 
k1      = fixpars(64); % transition A1 to A2 rate constant, 1/sec
km1     = fixpars(65); % transition A2 to A1 rate constant, 1/sec
alpha1  = fixpars(66); % Stretch sensing parameter for k1 and k?1, 1/um 
k2      = fixpars(67); % transition A2 to A3 rate constant, 1/sec
km2     = fixpars(68); % transition A3 to A2 rate constant, 1/sec
alpha2  = fixpars(69); % Stretch sensing parameter for k2 and k?2, 1/um  
K_D     = fixpars(70); % Used the values from Tewari etal JMCC (9/5 BM)
K_T     = fixpars(71); 
k3      = fixpars(72); % transition A3 to P rate constant, 1/sec
alpha3  = fixpars(73); % Stretch sensing parameter for k3, 1/um 
s3      = fixpars(74); % Strain at which k3 is minimum, um
k_on    = fixpars(75); % Campbell et al, Biophysical Journal 2018
k_off_LV   = fixpars(76);
k_off_RV   = fixpars(77);
K_coop  = fixpars(78); % Campbell et al, Biophysical Journal 2018,
ksr     = fixpars(79);
kmsr    = fixpars(80); % fit to the invitro data
kforce  = fixpars(81); 

%% Variables 

% Axial distance of midwall junction (cm)
xm_LV  = x(1); 
xm_SEP = x(2); 
xm_RV  = x(3);

% Radial distance of midwall junction (cm)
ym = x(4); 

% Contractile element length (um)
Lsc_LV  = x(5); 
Lsc_SEP = x(6); 
Lsc_RV  = x(7); 

% Volumes (cm^3) 
V_LV = x(8); 
V_RV = x(9);
V_SV = x(10);
V_PV = x(11); 
V_SA = x(12); 
V_PA = x(13); 

% Calcium moments 
P1_0_LV = x(14); % 0th moment state A1, LV
P1_1_LV = x(15); % 1st moment state A1, LV
P1_2_LV = x(16); % 2nd moment state A1, LV
P2_0_LV = x(17); % 0th moment state A2, LV
P2_1_LV = x(18); % 1st moment state A2, LV
P2_2_LV = x(19); % 2nd moment state A2, LV
P3_0_LV = x(20); % 0th moment state A3, LV
P3_1_LV = x(21); % 1st moment state A3, LV
P3_2_LV = x(22); % 2nd moment state A3, LV
N_LV = x(23);
U_NR_LV = x(24);
P1_0_SEP = x(25); % 0th moment state A1, LV
P1_1_SEP = x(26); % 1st moment state A1, LV
P1_2_SEP = x(27); % 2nd moment state A1, LV
P2_0_SEP = x(28); % 0th moment state A2, LV
P2_1_SEP= x(29);  % 1st moment state A2, LV
P2_2_SEP = x(30); % 2nd moment state A2, LV
P3_0_SEP = x(31); % 0th moment state A3, LV
P3_1_SEP = x(32); % 1st moment state A3, LV
P3_2_SEP = x(33); % 2nd moment state A3, LV
N_SEP = x(34);
U_NR_SEP = x(35);
P1_0_RV = x(36); % 0th moment state A1, LV
P1_1_RV = x(37); % 1st moment state A1, LV
P1_2_RV = x(38); % 2nd moment state A1, LV
P2_0_RV = x(39); % 0th moment state A2, LV
P2_1_RV = x(40); % 1st moment state A2, LV
P2_2_RV = x(41); % 2nd moment state A2, LV
P3_0_RV = x(42); % 0th moment state A3, LV
P3_1_RV = x(43); % 1st moment state A3, LV
P3_2_RV = x(44); % 2nd moment state A3, LV
N_RV = x(45);
U_NR_RV = x(46);

%% Correcting rate constants for metabolite levels in LV, SEP, and RV
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

%% Ca Calculations 

% Time to maximal systole 
TS_v = k_TS * T; 

% Time from maximal systole to relaxation 
TR_v = k_TR * T; 

Ca_amplitude = 2;

Ca_diastole = 0.1610;

tc = mod(t,T);
if tc >= 0 && tc < TS_v
    Ca_i = Ca_diastole + Ca_amplitude*0.5*(1 - cos(pi*(tc)/TS_v)); 
elseif tc >= TS_v && tc < TR_v + TS_v 
    Ca_i = Ca_diastole + Ca_amplitude*0.5*(1 + cos(pi*(tc - TS_v)/TR_v)); 
else
    Ca_i = Ca_diastole; 
end 
% figure(1111)
% hold on
% plot(t,Ca_i,'bx')

%% Heart model

% Volume of spherical cap formed by midwall surface (cm^3)
Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
Vm_RV  = (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 

% Surface area of midwall surface (cm^2) 
Am_LV  = pi * (xm_LV^2  + ym^2);
Am_SEP = pi * (xm_SEP^2 + ym^2); 
Am_RV  = pi * (xm_RV^2  + ym^2); 

% Curvature of midwall surface (cm^(-1))
Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
Cm_SEP =  2 * xm_SEP / (xm_SEP^2 + ym^2);
Cm_RV  =  2 * xm_RV  / (xm_RV^2  + ym^2);

% Ratio of wall thickness to midwall radius of curvature (dimensionless)
z_LV  = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV); 
z_SEP = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP); 
z_RV  = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV);

% Myofiber strain (dimensionless)
eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 

% Sarcomere length (cm)
Ls_LV  = Lsref * exp(eps_LV); 
Ls_SEP = Lsref * exp(eps_SEP); 
Ls_RV  = Lsref * exp(eps_RV); 

% % Collagen force
% sigma_collagen_LV  = PConcollagen*(exp(PExpcollagen*(Ls_LV - SLcollagen)) - 1).*(Ls_LV > SLcollagen);
% sigma_collagen_SEP = PConcollagen*(exp(PExpcollagen*(Ls_SEP - SLcollagen)) - 1).*(Ls_SEP > SLcollagen);
% sigma_collagen_RV  = PConcollagen*(exp(PExpcollagen*(Ls_RV - SLcollagen)) - 1).*(Ls_RV > SLcollagen);
% 
% % Passive Stress (kPa) 
% sigma_pas_LV_0  = 1*k_passive*(Ls_LV/2-L_rest_pas)  ;
% sigma_pas_SEP_0 = 1*k_passive*(Ls_SEP/2-L_rest_pas) ;
% sigma_pas_RV_0  = eta_sigma_pas*k_passive*(Ls_RV/2-L_rest_pas);

% sigma_pas_LV_old  = sigma_pas_LV_0  + sigma_collagen_LV ;
% sigma_pas_SEP_old = sigma_pas_SEP_0 + sigma_collagen_SEP;
% sigma_pas_RV_old  = sigma_pas_RV_0  + sigma_collagen_RV;

sigma_pas_LV  =  real(k_passive_LV * ((Ls_LV - Lsc0))^gamma); 
sigma_pas_SEP =  real(k_passive_LV * ((Ls_SEP - Lsc0))^gamma); 
sigma_pas_RV  =  real(k_passive_RV * ((Ls_RV - Lsc0))^gamma); 

% Sarcomere geometry (um) 
sovr_ze = min(L_thick*0.5, Lsc_LV*0.5);
sovr_cle = max(Lsc_LV*0.5 - (Lsc_LV-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_LV = L_sovr*2/(L_thick - L_hbare);
% sov_thin_LV = L_sovr/L_thin;

sovr_ze = min(L_thick*0.5, Lsc_SEP*0.5);
sovr_cle = max(Lsc_SEP*0.5 - (Lsc_SEP-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_SEP = L_sovr*2/(L_thick - L_hbare);
% sov_thin_SEP = L_sovr/L_thin;

sovr_ze = min(L_thick*0.5, Lsc_RV*0.5);
sovr_cle = max(Lsc_RV*0.5 - (Lsc_RV-L_thin),L_hbare*0.5);
L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
N_overlap_RV = L_sovr*2/(L_thick - L_hbare); % unitless 
% sov_thin_RV = L_sovr/L_thin;

% Active stress (kPa)
sigma_XB_LV  = N_overlap_LV*(kstiff2_LV*deltaR*(P3_0_LV) + kstiff1_LV*(P2_1_LV + P3_1_LV)); % mmHg * normalised force
sigma_XB_SEP = N_overlap_SEP*(kstiff2_LV*deltaR*(P3_0_SEP) + kstiff1_LV*(P2_1_SEP+P3_1_SEP)); % mmHg * normalised force
sigma_XB_RV  = N_overlap_RV*(kstiff2_RV*deltaR*(P3_0_RV) + kstiff1_RV*(P2_1_RV+P3_1_RV)); % mmHg * normalised force

% Total stress (kPa)
sigma_LV = -Kse*(Lsc_LV - Ls_LV);
sigma_SEP = -Kse*(Lsc_SEP - Ls_SEP);
sigma_RV = -Kse*(Lsc_RV - Ls_RV);

% Representative midwall tension (kPa cm)
Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

% Axial midwall tension component (kPa cm)
Tx_LV  = - Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 

% Radial midwall tension component (kPa cm)
Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2); 
Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2); 
Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);

% Ventricular pressure (kPa)
ptrans_LV = 2 * Tx_LV / ym; 
ptrans_RV = 2 * Tx_RV / ym; 
P_LV      = -ptrans_LV; 
P_RV      = ptrans_RV; 

% Calculate Pericardial Pressure 
Vh = V_LV + V_RV; 
P_peri = 0; % exp(s*(Vh/Vh0 - 1))/7.5; % convert mmHg to kPa 
P_LV = P_peri + P_LV; 
P_RV = P_peri + P_RV; 

%% Lumped circulatory model 
% Pressure (kPa)
P_SA = V_SA / C_SA;
P_SV = V_SV / C_SV; 
P_PA = V_PA / C_PA;
P_PV = V_PV / C_PV;  

% Flow 
% % When aortic valve is closed 
% Q_a = 0; 
% P_SA = (R_SA*V_SA + C_SA*P_SV*R_tSA + C_SA*Q_a*R_SA*R_tSA)/(C_SA*(R_SA + R_tSA)); 
% Q_SA = (V_SA - C_SA*P_SV + C_SA*Q_a*R_tSA)/(C_SA*(R_SA + R_tSA)); 
% 
% % When aortic valve is open 
% if (P_SA < P_LV) * (V_LV > 0)
%     Q_a = -(R_SA*V_SA - C_SA*P_LV*R_SA - C_SA*P_LV*R_tSA + C_SA*P_SV*R_tSA)/(C_SA*(R_a_vlv*R_SA + R_a_vlv*R_tSA + R_SA*R_tSA)); 
%     Q_SA = (R_a_vlv*V_SA - C_SA*P_SV*R_a_vlv + C_SA*P_LV*R_tSA - C_SA*P_SV*R_tSA)/(C_SA*(R_a_vlv*R_SA + R_a_vlv*R_tSA + R_SA*R_tSA)); 
%     P_SA = (R_a_vlv*R_SA*V_SA + C_SA*P_SV*R_a_vlv*R_tSA + C_SA*P_LV*R_SA*R_tSA)/(C_SA*(R_a_vlv*R_SA + R_a_vlv*R_tSA + R_SA*R_tSA)); 
% end
% Q_a = max(Q_a,0); 
% 
% % When pulmonary valve is closed 
% Q_p = 0; 
% P_PA = (R_PA*V_PA + C_PA*P_PV*R_tPA + C_PA*Q_p*R_PA*R_tPA)/(C_PA*(R_PA + R_tPA)); 
% Q_PA = (V_PA - C_PA*P_PV + C_PA*Q_p*R_tPA)/(C_PA*(R_PA + R_tPA)); 
% 
% % When pulmonary valve is open 
% if (P_PA < P_RV) * (V_RV > 0) 
%     Q_p = -(R_PA*V_PA - C_PA*P_RV*R_PA + C_PA*P_PV*R_tPA - C_PA*P_RV*R_tPA)/(C_PA*(R_PA*R_p_vlv + R_PA*R_tPA + R_p_vlv*R_tPA)); 
%     Q_PA = (R_p_vlv*V_PA - C_PA*P_PV*R_p_vlv - C_PA*P_PV*R_tPA + C_PA*P_RV*R_tPA)/(C_PA*(R_PA*R_p_vlv + R_PA*R_tPA + R_p_vlv*R_tPA)); 
%     P_PA = (R_PA*R_p_vlv*V_PA + C_PA*P_PV*R_p_vlv*R_tPA + C_PA*P_RV*R_PA*R_tPA)/(C_PA*(R_PA*R_p_vlv + R_PA*R_tPA + R_p_vlv*R_tPA));  
% end 
% Q_p = max(Q_p, 0); 

Q_m = max((P_PV - P_LV) / R_m_vlv, 0); 
Q_t = max((P_SV - P_RV) / R_t_vlv, 0); 
Q_a = max((P_LV - P_SA) / R_a_vlv, 0); 
Q_p = max((P_RV - P_PA) / R_p_vlv, 0); 

Q_SA = (P_SA - P_SV) / R_SA;
Q_PA = (P_PA - P_PV) / R_PA;

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
Joff_LV = k_off_LV*P0_LV*(1 + K_coop*N_LV);
% SEP 
P0_SEP = 1 - N_SEP - P1_0_SEP - P2_0_SEP - P3_0_SEP;  % DAB 10/8/2019

U_SR_SEP = 1 - U_NR_SEP;
Jon_SEP = k_on*Ca_i*N_SEP*(1 + K_coop*(1 - N_SEP));
Joff_SEP = k_off_LV*P0_SEP*(1 + K_coop*N_SEP);
% RV
P0_RV = 1.0 - N_RV - P1_0_RV - P2_0_RV - P3_0_RV;  % DAB 10/8/2019

U_SR_RV = 1.0 - U_NR_RV;
Jon_RV = k_on*Ca_i*N_RV*(1 + K_coop*(1 - N_RV));
Joff_RV = k_off_RV*P0_RV*(1 + K_coop*N_RV);

%% ODEs
% dxdt 
% 1 - 4 (triseg) 
dxm_LV  = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
dxm_SEP = Tx_LV + Tx_SEP + Tx_RV;
dxm_RV  = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
dym     = Ty_LV + Ty_SEP + Ty_RV; 

% 5 - 7 (Myofiber Mechanics)
dLsc_LV = (sigma_LV - sigma_pas_LV - sigma_XB_LV)/eta;
dLsc_SEP = (sigma_SEP - sigma_pas_SEP - sigma_XB_SEP)/eta;
dLsc_RV = (sigma_RV - sigma_pas_RV - sigma_XB_RV)/eta;

% 8 - 18 (xb, LV) 
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
dU_NR_LV = ksr * (1 + kforce * sigma_XB_LV) * U_SR_LV - kmsr*U_NR_LV ; % DAB 10/8

% 19 - 29 (xb, SEP) 
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
dU_NR_SEP = ksr*(1 + kforce * sigma_XB_SEP) * U_SR_SEP - kmsr*U_NR_SEP ; % DAB 10/8

% 30 - 40 (xb, RV) 
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
dU_NR_RV = ksr*(1 + kforce * sigma_XB_RV) * U_SR_RV - kmsr * U_NR_RV ; % DAB 10/8

% 41 - 46 (Lumped circulation model variables) 
dV_LV = Q_m - Q_a; 
dV_SA = Q_a - Q_SA; 
dV_SV = Q_SA - Q_t; 
dV_RV = Q_t - Q_p; 
dV_PA = Q_p - Q_PA; 
dV_PV = Q_PA - Q_m; 

dxdt = [dxm_LV; dxm_SEP; dxm_RV; dym;           % 1-4
    dLsc_LV; dLsc_SEP; dLsc_RV;                 % 5-7   
    dV_LV; dV_RV; dV_SV; dV_PV; dV_SA; dV_PA;   % 8-13
    dP1_0_LV; dP1_1_LV; dP1_2_LV;               % 14-16
    dP2_0_LV; dP2_1_LV; dP2_2_LV;               % 17-19
    dP3_0_LV; dP3_1_LV; dP3_2_LV;               % 20-22
    dN_LV; dU_NR_LV;                            % 23-24
    dP1_0_SEP; dP1_1_SEP; dP1_2_SEP;            % 25-27
    dP2_0_SEP; dP2_1_SEP; dP2_2_SEP;            % 28-30
    dP3_0_SEP; dP3_1_SEP; dP3_2_SEP;            % 31-33
    dN_SEP; dU_NR_SEP;                          % 34-35
    dP1_0_RV; dP1_1_RV; dP1_2_RV;               % 36-38
    dP2_0_RV; dP2_1_RV; dP2_2_RV;               % 39-41
    dP3_0_RV; dP3_1_RV; dP3_2_RV;               % 42-44
    dN_RV; dU_NR_RV;                            % 45-46
    ];  


outputs = [P_LV; P_SA; P_SV; P_RV; P_PA; P_PV;  % 1-6
    Vm_LV; Vm_SEP; Vm_RV;                       % 7-9
    Am_LV; Am_SEP; Am_RV;                       % 10-12
    Cm_LV; Cm_SEP; Cm_RV;                       % 13-15
    eps_LV; eps_SEP; eps_RV;                    % 16-18
    sigma_pas_LV; sigma_pas_SEP; sigma_pas_RV;  % 19-21
    sigma_XB_LV; sigma_XB_SEP; sigma_XB_RV;  % 22-24
    sigma_LV; sigma_SEP; sigma_RV;              % 25-27
    Q_m; Q_a; Q_t; Q_p;                         % 28-31
    Q_SA; Q_PA;                                 % 32-33
    Tm_LV; Tm_SEP; Tm_RV;                       % 34-36
    P_peri; Ca_i;                               % 37-38
    P0_LV; P0_SEP; P0_RV                        % 
    ];

end 