function f = triseg(x,pars,data,init)

%{ 
    This function contains the equations to calculate consistent initial
    conditions for the DAE. 
    Inputs: 
    x           - vector of states 
    pars        - vector of adjustable parameters 
    data        - input data structure with data and global parameters 
    init        - vector of initial conditions 
    Outputs: 
    f           - vector of equations for the root finder 
%} 

x = exp(x); 

HR = data.HR; 

fixpars = data.fixpars;

%% Adjustable Parameters
% Midwall area 
Amref_LV  = pars(11); 
Amref_SEP = pars(12); 
Amref_RV  = pars(13); 

% Midwall volume 
Vw_LV  = pars(14); 
Vw_SEP = pars(15); 
Vw_RV  = pars(16); 

%% Fixed parameters

% Sarcomere length parameters 
Lsref   = fixpars(13);
Kse    = fixpars(27); % series element elastance, mmHg/micron to kPa/cm

%% Variables 

xm_LV  = x(1); 
xm_SEP = x(2); 
xm_RV  = x(3); 
ym     = x(4); 

% Contractile element length 
Lsc_LV  = init(5); 
Lsc_SEP = init(6); 
Lsc_RV  = init(7); 

% Volumes 
V_LV = init(8); 
V_RV = init(9);

%% Heart and sarcomere model 

% Volume of spherical cap formed by midwall surface (m^3)
Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
Vm_SEP =  (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
Vm_RV  =  (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 

% Surface area of midwall surface (m^2) 
Am_LV  = pi * (xm_LV^2  + ym^2);
Am_SEP = pi * (xm_SEP^2 + ym^2); 
Am_RV  = pi * (xm_RV^2  + ym^2); 

% Curvature of midwall surface (m^(-1))
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

% Sarcomere length (m)
Ls_LV  = Lsref * exp(eps_LV); 
Ls_SEP = Lsref * exp(eps_SEP); 
Ls_RV  = Lsref * exp(eps_RV); 

% Total stress (kPa)
sigma_LV  = -Kse*(Lsc_LV - Ls_LV);
sigma_SEP = -Kse*(Lsc_SEP - Ls_SEP);
sigma_RV  = -Kse*(Lsc_RV - Ls_RV);

% Representative midwall tension (mmHg cm^(-2))
Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

% Axial midwall tension component 
Tx_LV  = - Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
Tx_SEP =   Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
Tx_RV  =   Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 

% Radial midwall tension component 
Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2); 
Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2); 
Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);

%% System of equations 

f1 = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
f2 = Tx_LV + Tx_SEP + Tx_RV;
f3 = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
f4 = Ty_LV + Ty_SEP + Ty_RV; 

% f5 = (Ls_LV  - Lsc_LV)  /Lse_iso - 1;
% f6 = (Ls_SEP - Lsc_SEP) /Lse_iso - 1;
% f7 = (Ls_RV  - Lsc_RV)  /Lse_iso - 1;

% f = [f1; f2; f3; f4; f5; f6; f7]; 
f = [f1; f2; f3; f4]; 

end 

