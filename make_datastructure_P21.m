function [data] = make_datastructure_P21(animal_id,table)

% Get the animal for the current animal
loc = find(table.AnimalID == animal_id);
t_current = table(loc,:);

T = t_current.HeartPeriod;

BW = t_current.Bodyweight;

WR_LVandSeptum = t_current.LVAndSeptumWeightRatio;
WR_RV = t_current.RVWeightRatio;

SV = t_current.StrokeVolume;

ESP_LV = t_current.LVESP;
EDP_LV = t_current.LVEDP;

ESP_RV = t_current.RVESP;
EDP_RV = t_current.RVEDP;

ESV_LV = t_current.LVESV;
ESV_RV = t_current.RVESV;

%% Calculate Data
% End-Diastolic Volume 
EDV_LV = ESV_LV + SV;
EDV_RV = ESV_RV + SV; 

% Set tspan to solve over the first 20 beats 
dt    = 0.001; 
tspan = 0:dt:20*T;

% Heart rate (beats/s --> beats/min)  
HR = 1/T*60;

% Total blood volume (mL = cm^3)
Vtot = BW*60*1e-3; % 50g bw and 60mL/kg of bw

% Cardiac Output (uL/min to cm^3 s^(-1))
CO_LV =  (HR*SV)/60/1000;

% Blood pressure in mmHg 
SPbar = ESP_LV/1.05; % ESP_LV + .5; 
DPbar = SPbar - 40; 

% Ventricular wall weights (mili g --> g) 
Weight_LV_and_SEP = WR_LVandSeptum*BW*1e-3;
Weight_RV         = WR_RV*BW*1e-3; 

% Ventricular wall volumes = = weight (g) / 1.055(g/mL) 
Wall_Volume_LV_and_SEP = Weight_LV_and_SEP/1.055; 
Wall_Volume_RV         = Weight_RV/1.055; 

%% Save in a data structure 
data.SPbar  = SPbar / 7.5; 
data.DPbar  = DPbar / 7.5; 
data.HR     = HR; 
data.T      = T; 
data.Vtot   = Vtot; 
data.CO     = CO_LV ; 
data.tspan  = tspan; 
data.dt     = dt;
data.EDP_LV = EDP_LV / 7.5; 
data.EDV_LV = EDV_LV * 1e-3;
data.ESP_LV = ESP_LV / 7.5; 
data.ESV_LV = ESV_LV * 1e-3; 
data.EDP_RV = EDP_RV / 7.5; 
data.EDV_RV = EDV_RV * 1e-3; 
data.ESP_RV = ESP_RV / 7.5; 
data.ESV_RV = ESV_RV * 1e-3; 

data.Wall_Volume_LV_and_SEP = Wall_Volume_LV_and_SEP;
data.Wall_Volume_RV         = Wall_Volume_RV; 

data.HR = HR; 
end