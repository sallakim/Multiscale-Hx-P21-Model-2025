% Get Table 5
clear 
clc

addpath opt_pars model outputs

nx_animal_ids = [11 12 51 52 54 55 56];
hx_animal_ids = [1 4 5 7 8 9 10 57 58 59 61 62]; 

counter = 0 ;
%%
for i = 1:length(nx_animal_ids)
animal_id = nx_animal_ids(i);
filename = sprintf('outputs_Nx%d.mat', animal_id);
filename2 = sprintf('opt_pars_Nx%d.mat', animal_id);
load(filename)
load(filename2,'pars_opt')
optpars = exp(pars_opt);

V_LV = outputs.volumes.V_LV; % mL
V_RV = outputs.volumes.V_RV;
P_LV = outputs.pressures.P_LV; % mmHg 
P_RV = outputs.pressures.P_RV;

[EDV, EDP, ESV, ESP] = getEDESvals(V_LV,V_RV,P_LV,P_RV); % LV, RV 

EDV_LV_mat(i) = EDV(1); EDV_RV_mat(i) = EDV(2); 

ESV_LV_mat(i) = ESV(1); ESV_RV_mat(i) = ESV(2); 

EDP_LV_mat(i) = EDP(1); EDP_RV_mat(i) = EDP(2); 

ESP_LV_mat(i) = ESP(1); ESP_RV_mat(i) = ESP(2); 

% Get data 

table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
data = make_datastructure_P21(animal_id,table);

EDV_LV_data = data.EDV_LV; EDV_RV_data = data.EDV_RV; % mL
ESV_LV_data = data.ESV_LV; ESV_RV_data = data.ESV_RV; 
EDP_LV_data = data.EDP_LV * 7.5; EDP_RV_data = data.EDP_RV * 7.5; % 
ESP_LV_data = data.ESP_LV * 7.5; ESP_RV_data = data.ESP_RV * 7.5;

EDV_LV_data_mat(i) = EDV_LV_data; EDV_RV_data_mat(i) = EDV_RV_data; 

ESV_LV_data_mat(i) = ESV_LV_data; ESV_RV_data_mat(i) = ESV_RV_data; 

EDP_LV_data_mat(i) = EDP_LV_data; EDP_RV_data_mat(i) = EDP_RV_data; 

ESP_LV_data_mat(i) = ESP_LV_data; ESP_RV_data_mat(i) = ESP_RV_data; 

R_SA(i) = optpars(5); R_PA(i) = optpars(6); 

counter = i; 
end

% For EDV_LV 
% numerator = sum( (EDV_LV_data_mat - EDV_LV_mat).^2 );
% denominator = sum( (EDV_LV_data_mat - mean(EDV_LV_data_mat)).^2 );
% R_squared_EDV_LV = 1 - numerator/denominator
[RHO_S,PVAL_S] = corr(EDV_LV_data_mat',EDV_LV_mat','Type','Spearman');
fprintf(['EDV_LV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDV_LV_mat), std(EDV_LV_mat),RHO_S)

% For ESV_LV
% numerator_ESV_LV = sum( (ESV_LV_data_mat - ESV_LV_mat).^2 );
% denominator_ESV_LV = sum( (ESV_LV_data_mat - mean(ESV_LV_data_mat)).^2 );
% R_squared_ESV_LV = 1 - numerator_ESV_LV/denominator_ESV_LV
[RHO_S,PVAL_S] = corr(ESV_LV_data_mat',ESV_LV_mat','Type','Spearman');
fprintf(['ESV_LV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(ESV_LV_mat), std(ESV_LV_mat),RHO_S)

% For EDV_RV
% numerator_EDV_RV = sum( (EDV_RV_data_mat - EDV_RV_mat).^2 );
% denominator_EDV_RV = sum( (EDV_RV_data_mat - mean(EDV_RV_data_mat)).^2 );
% R_squared_EDV_RV = 1 - numerator_EDV_RV/denominator_EDV_RV
[RHO_S,PVAL_S] = corr(EDV_RV_data_mat',EDV_RV_mat','Type','Spearman');
fprintf(['EDV_RV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDV_RV_mat), std(EDV_RV_mat),RHO_S)

% For ESV_RV
% numerator_ESV_RV = sum( (ESV_RV_data_mat - ESV_RV_mat).^2 );
% denominator_ESV_RV = sum( (ESV_RV_data_mat - mean(ESV_RV_data_mat)).^2 );
% R_squared_ESV_RV = 1 - numerator_ESV_RV/denominator_ESV_RV
[RHO_S,PVAL_S] = corr(ESV_RV_data_mat',ESV_RV_mat','Type','Spearman');
fprintf(['ESV_RV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(ESV_RV_mat), std(ESV_RV_mat),RHO_S)

% For EDP_LV
% numerator_EDP_LV = sum( (EDP_LV_data_mat - EDP_LV_mat).^2 );
% denominator_EDP_LV = sum( (EDP_LV_data_mat - mean(EDP_LV_data_mat)).^2 );
% R_squared_EDP_LV = 1 - numerator_EDP_LV/denominator_EDP_LV
[RHO_S,PVAL_S] = corr(EDP_LV_data_mat',EDP_LV_mat','Type','Spearman');
fprintf(['EDP_LV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDP_LV_mat), std(EDP_LV_mat),RHO_S)

% For ESP_LV
% numerator_ESP_LV = sum( (ESP_LV_data_mat - ESP_LV_mat).^2 );
% denominator_ESP_LV = sum( (ESP_LV_data_mat - mean(ESP_LV_data_mat)).^2 );
% R_squared_ESP_LV = 1 - numerator_ESP_LV/denominator_ESP_LV
[RHO_S,PVAL_S] = corr(ESP_LV_data_mat',ESP_LV_mat','Type','Spearman');
fprintf(['ESP_LV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(ESP_LV_mat), std(ESP_LV_mat),RHO_S)

% For EDP_RV
% numerator_EDP_RV = sum( (EDP_RV_data_mat - EDP_RV_mat).^2 );
% denominator_EDP_RV = sum( (EDP_RV_data_mat - mean(EDP_RV_data_mat)).^2 );
% R_squared_EDP_RV = 1 - numerator_EDP_RV/denominator_EDP_RV
[RHO_S,PVAL_S] = corr(EDP_RV_data_mat',EDP_RV_mat','Type','Spearman');
fprintf(['EDP_RV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDP_RV_mat), std(EDP_RV_mat),RHO_S)

% For ESP_RV
% numerator_ESP_RV = sum( (ESP_RV_data_mat - ESP_RV_mat).^2 );
% denominator_ESP_RV = sum( (ESP_RV_data_mat - mean(ESP_RV_data_mat)).^2 );
% R_squared_ESP_RV = 1 - numerator_ESP_RV/denominator_ESP_RV
[RHO_S,PVAL_S] = corr(ESP_RV_data_mat',ESP_RV_mat','Type','Spearman');
fprintf(['ESP_RV_nx, Mean model: %d, std: %d, R2: %d \n' ...
    '\n'],mean(ESP_RV_mat), std(ESP_RV_mat),RHO_S)


%%

for i = 1:length(hx_animal_ids)
animal_id = hx_animal_ids(i);
filename = sprintf('outputs_Hx%d.mat', animal_id);
filename2 = sprintf('opt_pars_Hx%d.mat', animal_id);
load(filename)
load(filename2,'pars_opt')
optpars = exp(pars_opt);

counter = counter + 1; 

V_LV = outputs.volumes.V_LV; % mL
V_RV = outputs.volumes.V_RV;
P_LV = outputs.pressures.P_LV; % mmHg 
P_RV = outputs.pressures.P_RV;

[EDV, EDP, ESV, ESP] = getEDESvals(V_LV,V_RV,P_LV,P_RV); % LV, RV 

EDV_LV_mat(counter) = EDV(1); EDV_RV_mat(counter) = EDV(2); 

ESV_LV_mat(counter) = ESV(1); ESV_RV_mat(counter) = ESV(2); 

EDP_LV_mat(counter) = EDP(1); EDP_RV_mat(counter) = EDP(2); 

ESP_LV_mat(counter) = ESP(1); ESP_RV_mat(counter) = ESP(2); 

% Get data 
EDV_LV_data = data.EDV_LV; EDV_RV_data = data.EDV_RV; % mL
ESV_LV_data = data.ESV_LV; ESV_RV_data = data.ESV_RV; 
EDP_LV_data = data.EDP_LV * 7.5; EDP_RV_data = data.EDP_RV * 7.5; % 
ESP_LV_data = data.ESP_LV * 7.5; ESP_RV_data = data.ESP_RV * 7.5;

EDV_LV_data_mat(counter) = EDV_LV_data; EDV_RV_data_mat(counter) = EDV_RV_data; 

ESV_LV_data_mat(counter) = ESV_LV_data; ESV_RV_data_mat(counter) = ESV_RV_data; 

EDP_LV_data_mat(counter) = EDP_LV_data; EDP_RV_data_mat(counter) = EDP_RV_data; 

ESP_LV_data_mat(counter) = ESP_LV_data; ESP_RV_data_mat(counter) = ESP_RV_data; 

end

% For EDV_LV 
% numerator = sum( (EDV_LV_data_mat - EDV_LV_mat).^2 );
% denominator = sum( (EDV_LV_data_mat - mean(EDV_LV_data_mat)).^2 );
% R_squared_EDV_LV = 1 - numerator/denominator
[RHO_S,PVAL_S] = corr(EDV_LV_data_mat',EDV_LV_mat','Type','Spearman');
fprintf(['EDV_LV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDV_LV_mat), std(EDV_LV_mat),RHO_S)

% For ESV_LV
% numerator_ESV_LV = sum( (ESV_LV_data_mat - ESV_LV_mat).^2 );
% denominator_ESV_LV = sum( (ESV_LV_data_mat - mean(ESV_LV_data_mat)).^2 );
% R_squared_ESV_LV = 1 - numerator_ESV_LV/denominator_ESV_LV
[RHO_S,PVAL_S] = corr(ESV_LV_data_mat',ESV_LV_mat','Type','Spearman');
fprintf(['ESV_LV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(ESV_LV_mat), std(ESV_LV_mat),RHO_S)

% For EDV_RV
% numerator_EDV_RV = sum( (EDV_RV_data_mat - EDV_RV_mat).^2 );
% denominator_EDV_RV = sum( (EDV_RV_data_mat - mean(EDV_RV_data_mat)).^2 );
% R_squared_EDV_RV = 1 - numerator_EDV_RV/denominator_EDV_RV
[RHO_S,PVAL_S] = corr(EDV_RV_data_mat',EDV_RV_mat','Type','Spearman');
fprintf(['EDV_RV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDV_RV_mat), std(EDV_RV_mat),RHO_S)

% For ESV_RV
% numerator_ESV_RV = sum( (ESV_RV_data_mat - ESV_RV_mat).^2 );
% denominator_ESV_RV = sum( (ESV_RV_data_mat - mean(ESV_RV_data_mat)).^2 );
% R_squared_ESV_RV = 1 - numerator_ESV_RV/denominator_ESV_RV
[RHO_S,PVAL_S] = corr(ESV_RV_data_mat',ESV_RV_mat','Type','Spearman');
fprintf(['ESV_RV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(ESV_RV_mat), std(ESV_RV_mat),RHO_S)

% For EDP_LV
% numerator_EDP_LV = sum( (EDP_LV_data_mat - EDP_LV_mat).^2 );
% denominator_EDP_LV = sum( (EDP_LV_data_mat - mean(EDP_LV_data_mat)).^2 );
% R_squared_EDP_LV = 1 - numerator_EDP_LV/denominator_EDP_LV
[RHO_S,PVAL_S] = corr(EDP_LV_data_mat',EDP_LV_mat','Type','Spearman');
fprintf(['EDP_LV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDP_LV_mat), std(EDP_LV_mat),RHO_S)

% For ESP_LV
% numerator_ESP_LV = sum( (ESP_LV_data_mat - ESP_LV_mat).^2 );
% denominator_ESP_LV = sum( (ESP_LV_data_mat - mean(ESP_LV_data_mat)).^2 );
% R_squared_ESP_LV = 1 - numerator_ESP_LV/denominator_ESP_LV
[RHO_S,PVAL_S] = corr(ESP_LV_data_mat',ESP_LV_mat','Type','Spearman');
fprintf(['ESP_LV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(ESP_LV_mat), std(ESP_LV_mat),RHO_S)

% For EDP_RV
% numerator_EDP_RV = sum( (EDP_RV_data_mat - EDP_RV_mat).^2 );
% denominator_EDP_RV = sum( (EDP_RV_data_mat - mean(EDP_RV_data_mat)).^2 );
% R_squared_EDP_RV = 1 - numerator_EDP_RV/denominator_EDP_RV
[RHO_S,PVAL_S] = corr(EDP_RV_data_mat',EDP_RV_mat','Type','Spearman');
fprintf(['EDP_RV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(EDP_RV_mat), std(EDP_RV_mat),RHO_S)

% For ESP_RV
% numerator_ESP_RV = sum( (ESP_RV_data_mat - ESP_RV_mat).^2 );
% denominator_ESP_RV = sum( (ESP_RV_data_mat - mean(ESP_RV_data_mat)).^2 );
% R_squared_ESP_RV = 1 - numerator_ESP_RV/denominator_ESP_RV
[RHO_S,PVAL_S] = corr(ESP_RV_data_mat',ESP_RV_mat','Type','Spearman');
fprintf(['ESP_RV_hx, Mean model: %d, std: %d, R2: %d \n' ...
    ],mean(ESP_RV_mat), std(ESP_RV_mat),RHO_S)
