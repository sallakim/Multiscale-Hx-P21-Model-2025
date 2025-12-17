% LSA Model Wrap 

function y = LSA_wrap(pars,data) 

outputs = model_sol(pars,data);

% Volume and Pressure timecourse data available 
V_LV = outputs.volumes.V_LV;
V_RV = outputs.volumes.V_RV; 
P_LV = outputs.pressures.P_LV; 
P_RV = outputs.pressures.P_RV; 

% [EDV, EDP, ESV, ESP] = getEDESvals(V_LV,V_RV,P_LV,P_RV); 
% 
% EDV_LV = EDV(1); EDV_RV = EDV(2);
% EDP_LV = EDP(1); EDP_RV = EDP(2);
% ESV_LV = ESV(1); ESV_RV = ESV(2);
% ESP_LV = ESP(1); ESP_RV = ESP(2);
% 
% SV_LV = EDV_LV - ESV_LV; 
% SV_RV = EDV_RV - ESV_RV; 
% 
% EF_LV = SV_LV/EDV_LV; 
% EF_RV = SV_RV/EDV_RV; 
% 
% y = [EDV_LV; ESV_LV; EDP_LV; ESP_LV;
%     EDV_RV; ESV_RV; EDP_RV; ESP_RV;
%     EF_LV; EF_RV]; 

y = [V_LV', V_RV', P_LV, P_RV];

end