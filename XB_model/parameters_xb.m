function [init,pars] = parameters_xb(data)

%% data

if isfield(data,'eta_k_on') == 1
    eta_k_on = data.eta_k_on; 
else
    eta_k_on = 1; 
end

if isfield(data,'eta_k_off') == 1
    eta_k_off = data.eta_k_off; 
else
    eta_k_off = 1; 
end

if isfield(data,'eta_kstiff') == 1
    eta_kstiff = data.eta_kstiff; 
else
    eta_kstiff = 1; 
end

if isfield(data,'eta_K_coop') == 1
    eta_K_coop = data.eta_K_coop; 
else
    eta_K_coop = 1; 
end

if isfield(data,'eta_k3') == 1
    eta_k3 = data.eta_k3; 
else
    eta_k3 = 1; 
end

if isfield(data,'eta_K_D') == 1
    eta_K_D = data.eta_K_D; 
else
    eta_K_D = 1; 
end

if isfield(data,'eta_K_T') == 1
    eta_K_T = data.eta_K_T; 
else
    eta_K_T = 1; 
end

if isfield(data,'eta_K_Pi') == 1
    eta_K_Pi = data.eta_K_Pi; 
else
    eta_K_Pi = 1; 
end

if isfield(data,'eta_k2') == 1
    eta_k2 = data.eta_k2; 
else
    eta_k2 = 1; 
end

if isfield(data,'eta_k_passive') == 1
    eta_k_passive = data.eta_k_passive; 
else
    eta_k_passive = 1; 
end

if isfield(data,'eta_kd') == 1
    eta_kd = data.eta_kd; 
else
    eta_kd = 1; 
end

if isfield(data,'eta_ka') == 1
    eta_ka = data.eta_ka; 
else
    eta_ka = 1; 
end

if isfield(data,'eta_k1') == 1
    eta_k1 = data.eta_k1; 
else
    eta_k1 = 1; 
end

if isfield(data,'eta_km1') == 1
    eta_km1 = data.eta_km1; 
else
    eta_km1 = 1; 
end

if isfield(data,'eta_km2') == 1
    eta_km2 = data.eta_km2; 
else
    eta_km2 = 1; 
end

if isfield(data,'eta_kforce') == 1
    eta_kforce = data.eta_kforce; 
else
    eta_kforce = 1; 
end

if isfield(data,'eta_ksr') == 1
    eta_ksr = data.eta_ksr; 
else
    eta_ksr = 1; 
end

if isfield(data,'eta_kmsr') == 1
    eta_kmsr = data.eta_kmsr; 
else
    eta_kmsr = 1; 
end

if isfield(data,'eta_alpha1') == 1
    eta_alpha1 = data.eta_alpha1; 
else
    eta_alpha1 = 1; 
end

if isfield(data,'eta_alpha2') == 1
    eta_alpha2 = data.eta_alpha2; 
else
    eta_alpha2 = 1; 
end

if isfield(data,'eta_alpha3') == 1
    eta_alpha3 = data.eta_alpha3; 
else
    eta_alpha3 = 1; 
end

if isfield(data,'eta_s3') == 1
    eta_s3 = data.eta_s3; 
else
    eta_s3 = 1; 
end

%% Sarcomere lengths (um)  
Lsref = 2.0; 
Lsc_LV_0  = Lsref;
Lsc_SEP_0 = Lsref;
Lsc_RV_0  = Lsref;

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
kstiff1 = (1.4 * 1.7351e+03); % kPa/um (9/5 BM)
kstiff2 = eta_kstiff*(1.4 * 4.5545e+04);      % kPa/um (9/5 BM)
k_passive = eta_k_passive*50 /7.5;          % kPa / um % for mean SHAM rat and TAC rat 1

alpha3  = eta_alpha3*0.1*59.3; % Stretch sensing parameter for k3, 1/um 

k_off = eta_k_off*1*723.8520;     % manually tuned parameter!

% transitions between super relaxed state and non relaxed state
ksr = eta_ksr*15;
kforce = eta_kforce*1.169;         % 1/kpa
kmsr =    eta_kmsr*50.032;       % fit to the invitro data

ka      = eta_ka*0.8*559.5568; % myosin-actin attach rate constant, 1/sec
kd      = eta_kd*304.6708; % myosin-actin detach rate constant, 1/sec
k1      = eta_k1*112.3727; % transition A1 to A2 rate constant, 1/sec

k3      = eta_k3*144.5586; % transition A3 to P rate constant, 1/sec

% Collagen Froce 
SLcollagen = 2.25;   % threshold for collagen activation,(um)
PConcollagen = 0.01; % contriubtion of collagen (unitless)
PExpcollagen = 70;   % expresion of collagen (unitless)
% Passive Stress
eta       = 1/7.5; % visconsity, kPa s /micron (convert from mmHg) 
L_rest_pas = 0.9;  % (um)  

% Sarcomere Geometry Parameters (um)
L_thick = 1.67; % Length of thick filament
L_hbare = 0.10; % Length of bare region of thick filament
L_thin  = 1.20; % Length of thin filament
deltaR  = 0.010; 

Kse    = 50000 / 7.5; % series element elastance, kPa/micron (convert from mmHg/micron)

% para 4
km1     = eta_km1*21.296;   % transition A2 to A1 rate constant, 1/sec
k2      = eta_k2*811.72;   % transition A2 to A3 rate constant, 1/sec
km2     = eta_km2*43.25;    % transition A3 to A2 rate constant, 1/sec
alpha1  = eta_alpha1*10.0;     % Stretch sensing parameter for k1 and k?1, 1/um 
alpha2  = eta_alpha2*9.1;      % Stretch sensing parameter for k2 and k?2, 1/um  
s3      = eta_s3*9.9e-3;   % Strain at which k3 is minimum, um

K_coop = eta_K_coop*9.6846;               % Campbell et al, Biophysical Journal 2018, unitless
k_on =   eta_k_on*101.1850;    % Campbell et al, Biophysical Journal 2018, 1/(uM sec) 

K_Pi = eta_K_Pi*4.00e3;  % mM to um
K_T = eta_K_T*0.4897; % 
K_D = eta_K_D*0.194;  % Used the values from Tewari etal JMCC (9/5 BM)

%% Outputs 
init = [
    Lsc_LV_0; Lsc_SEP_0; Lsc_RV_0;                  % 1-3
    P1_0_LV; P1_1_LV; P1_2_LV ;P2_0_LV; P2_1_LV; P2_2_LV; P3_0_LV;P3_1_LV; P3_2_LV; N_LV; U_NR_LV;   % 4-14
    P1_0_SEP;P1_1_SEP;P1_2_SEP;P2_0_SEP;P2_1_SEP;P2_2_SEP;P3_0_SEP;P3_1_SEP;P3_2_SEP;N_SEP;U_NR_SEP; % 15-25
    P1_0_RV; P1_1_RV; P1_2_RV; P2_0_RV; P2_1_RV; P2_2_RV; P3_0_RV; P3_1_RV; P3_2_RV; N_RV; U_NR_RV]; % 26-36

pars = [          
        SLcollagen; PConcollagen; PExpcollagen;         % 1-3
        eta; L_rest_pas;                                % 4-5
        kstiff1; kstiff2; k_passive;                    % 6-8
        L_thick; L_hbare; L_thin; deltaR;               % 9-12
        alpha1; alpha2; alpha3;                         % 13-15
        s3;                                             % 16
        km1;  km2;                                      % 17-18
        k1; k2; k3;                                     % 19-21
        ka; kd;                                         % 22-23
        ksr; kforce; kmsr;                              % 24-26                           
        k_on;                                           % 27
        K_coop; Kse;                                    % 28-29
        K_Pi; K_T; K_D;                                 % 30-32
        k_off;                                          % 33
        ]; 
