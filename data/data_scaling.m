%% 
clear all

addpath raw_data\smoothed_data\

condition_flags = {'Hx f'};
load Hx59_LV_reshape.mat 
V_LV_avg = V_ven_avg; V_LV_stack = V_ven_stack; P_LV_avg = P_ven_avg; P_LV_stack = P_ven_stack; T_LV = T;
load Hx59_RV_reshape.mat
V_RV_avg = V_ven_avg; V_RV_stack = V_ven_stack; P_RV_avg = P_ven_avg; P_RV_stack = P_ven_stack; T_RV = T;
filename = 'Hx_delete.mat';

%% Data 
V_LV_avg_min = min(V_LV_avg); V_RV_avg_min = min(V_RV_avg); 
V_LV_avg_max = max(V_LV_avg); V_RV_avg_max = max(V_RV_avg);
P_LV_avg_min = min(P_LV_avg); P_RV_avg_min = min(P_RV_avg); 
T = (T_LV+T_RV)/2; T_span = linspace(0,T,50);
P_LV_mean = P_LV_avg; P_LV_mult = P_LV_stack; 
P_RV_mean = P_RV_avg; P_RV_mult = P_RV_stack; 
SV_common = max(V_LV_avg)-min(V_LV_avg)

%% Normalize
% LV
V_LV_norm = ( V_LV_avg - V_LV_avg_min )/( V_LV_avg_max - V_LV_avg_min );
V_LV_stack_norm = ( V_LV_stack - V_LV_avg_min )./( V_LV_avg_max - V_LV_avg_min );
% RV 
V_RV_norm = ( V_RV_avg - V_RV_avg_min )/( V_RV_avg_max - V_RV_avg_min );
V_RV_stack_norm = ( V_RV_stack - V_RV_avg_min )./( V_RV_avg_max - V_RV_avg_min );

%% Scaled to common stroke and shift back 
% LV 
V_LV_mean = V_LV_norm*SV_common + V_LV_avg_min; 
V_LV_mult = V_LV_stack_norm*SV_common + V_LV_avg_min; 
% RV 
V_RV_mean = V_RV_norm*SV_common + V_RV_avg_min; 
V_RV_mult = V_RV_stack_norm*SV_common + V_RV_avg_min; 

SV = SV_common; 

%% ESV adjustment 
if min(V_RV_mean) < 0 && min(V_LV_mean)> 0
    ESV_add = min(V_LV_mean)/3;
    V_RV_mean =  V_RV_mean + (0 - min(V_RV_mean))+ min(V_LV_mean) + ESV_add; 
elseif min(V_LV_mean) < 0 && min(V_RV_mean)> 0
    ESV_add = 0.7*min(V_RV_mean);
    V_LV_mean = V_LV_mean + (0-min(V_LV_mean)) + ESV_add;
elseif min(V_LV_mean) < 0 && min(V_RV_mean)< 0
    print('two negative volumes')
end

%% Pressure adjustment 
for i_exp = 1
condition_flag = condition_flags{i_exp};
switch condition_flag
    case 'Nx m' 
        if P_LV_avg_min < 0 
            P_LV_avg_min_new = 3.15-1; P_LV_mult = P_LV_mult + P_LV_avg_min_new + abs(P_LV_avg_min); P_LV_mean = P_LV_mean + P_LV_avg_min_new + abs(P_LV_avg_min);
        end
        if P_RV_avg_min < 0
            P_RV_avg_min_new = 2.55-1; P_RV_mult = P_RV_mult + P_RV_avg_min_new + abs(P_RV_avg_min); P_RV_mean = P_RV_mean + P_RV_avg_min_new + abs(P_RV_avg_min);
        end
    case 'Nx f'
        if P_LV_avg_min < 0 
            P_LV_avg_min_new = 3.35-1; P_LV_mult = P_LV_mult + P_LV_avg_min_new + abs(P_LV_avg_min); P_LV_mean = P_LV_mean + P_LV_avg_min_new + abs(P_LV_avg_min);
        end
        if P_RV_avg_min < 0
            P_RV_avg_min_new = 1.30-1; P_RV_mult = P_RV_mult + P_RV_avg_min_new + abs(P_RV_avg_min); P_RV_mean = P_RV_mean + P_RV_avg_min_new + abs(P_RV_avg_min);
        end
    case 'Hx m' 
        if P_LV_avg_min < 0 
            P_LV_avg_min_new = 4.8-1; P_LV_mult = P_LV_mult + P_LV_avg_min_new + abs(P_LV_avg_min); P_LV_mean = P_LV_mean + P_LV_avg_min_new + abs(P_LV_avg_min);
        end
        if P_RV_avg_min < 0
            P_RV_avg_min_new = 2.70-1; P_RV_mult = P_RV_mult + P_RV_avg_min_new + abs(P_RV_avg_min); P_RV_mean = P_RV_mean + P_RV_avg_min_new + abs(P_RV_avg_min);
        end
    case 'Hx f'
        if P_LV_avg_min < 0 
            P_LV_avg_min_new = 3.30-1; P_LV_mult = P_LV_mult + P_LV_avg_min_new + abs(P_LV_avg_min); P_LV_mean = P_LV_mean + P_LV_avg_min_new + abs(P_LV_avg_min);
        end
        if P_RV_avg_min < 0
            P_RV_avg_min_new = 3.40-1; P_RV_mult = P_RV_mult + P_RV_avg_min_new + abs(P_RV_avg_min); P_RV_mean = P_RV_mean + P_RV_avg_min_new + abs(P_RV_avg_min);
        end
end
end

%% Save 
save(filename,'T','SV','V_LV_mean','V_LV_mult','V_RV_mean','V_RV_mult','P_LV_mean','P_LV_mult','P_RV_mean','P_RV_mult')

minmax_LV = [min(min(V_LV_mean));max(max(V_LV_mean))];
minmax_RV = [min(min(V_RV_mean));max(max(V_RV_mean))];
A = [V_LV_mean;V_RV_mean;minmax_LV;minmax_RV;V_LV_avg;V_RV_avg];
minx = round(min(A)-0.001,3);
maxx = round(max(A)+0.001,3);

figure(1)
clf
subplot(1,2,1);
hold on
title(strcat(condition_flag),'LV') 
plot(T_span,V_LV_mean)
plot(T_span,V_LV_avg,'--')
ylim([minx,maxx])

subplot(1,2,2); 
hold on
title(strcat(condition_flag),'RV')
plot(T_span,V_RV_mean)
plot(T_span,V_RV_avg,'--')
ylim([minx,maxx])

figure(2)
clf
subplot(2,1,1);
hold on
title(strcat(condition_flag),'LV') 
plot(V_LV_mean,P_LV_mean)
plot(V_LV_avg,P_LV_avg,'--')
% ylim([0,round((max(P_LV_mean)+5))])
xlim([minx maxx])

subplot(2,1,2); 
hold on
title(strcat(condition_flag),'RV')
plot(V_RV_mean,P_RV_mean)
plot(V_RV_avg,P_RV_avg,'--')
% ylim([0,round((max(P_LV_mean)+5))])
xlim([minx maxx])

figure(3)
clf
subplot(1,2,1);
hold on
title(strcat(condition_flag),'LV') 
plot(T_span,P_LV_mean)
plot(T_span,P_LV_avg,'--')

subplot(1,2,2); 
hold on
title(strcat(condition_flag),'RV')
plot(T_span,P_RV_mean)
plot(T_span,P_RV_avg,'--')
