% Construct an "average" signal from the data
clear; clc; close all;
%%
table = readtable('Hx62_RV_channel_1_2.txt');
filename = 'Hx62_RV_smooth.mat';
data = table2array(table);
% table1 = readtable('Rat58_RV_channel1.txt');
% table2 = readtable('Rat58_RV_channel2.txt');
% data = table2array(table1);
% data(:,3:4) = table2array(table2);
nt = 50;

t_vals = data(:,1);
p_ventricle = data(:,2);   
V_ventricle = data(:,3);
dt = (t_vals(2)-t_vals(1))./1000; %ms --> s


smootherfac = 0.05;%0.05;
[p_ven_sm,pven_wind_size] = smoothdata(p_ventricle,'gaussian','SmoothingFactor',smootherfac);
[V_ven_sm,vven_wind_size] = smoothdata(V_ventricle,'gaussian','SmoothingFactor',smootherfac);
% [p_SA,psa_wind_size] = smoothdata(p_SA,'gaussian','SmoothingFactor',smootherfac);


% LV
% Use the full volume signal
%     [pks,loc] = findpeaks(V_LV,'MinPeakDistance',90);
% Or use the second derivative in pressure
temp = diff(diff(p_ven_sm)); 
[pks,loc] = findpeaks(temp,'MinPeakDistance',150);
loc = loc-1; % WAS minus 2


V_ven_stack = zeros(nt,length(pks)-1);
P_ven_stack = zeros(nt,length(pks)-1);
% PSA_stack = zeros(nt,length(pks)-1);

for j=1:length(loc)-1
    % Fine the signal points, then interpolate to be the same length for
    % averaging purposes
    temp = V_ventricle(loc(j):loc(j+1));
    V_ven_stack(:,j) = interp1(temp,linspace(1,length(temp),nt));

    temp = p_ventricle(loc(j):loc(j+1));   
    P_ven_stack(:,j) = interp1(temp,linspace(1,length(temp),nt));

%     temp = p_SA(loc(j):loc(j+1));
%     PSA_stack(:,j) = interp1(temp,linspace(1,length(temp),nt));

end
T_all = diff(t_vals(loc)/1000);
T = mean(T_all);


V_ven_avg = mean(V_ven_stack,2);
P_ven_avg = mean(P_ven_stack,2);
% P_SA_avg = mean(PSA_stack,2);

save (filename)

figure(1);clf;
subplot(1,2,1); hold on;
plot(V_ven_stack); plot(V_ven_avg,'k','LineWidth',3);
subplot(1,2,2); hold on;
plot(P_ven_stack); plot(P_ven_avg,'k','LineWidth',3);
set(gca,'FontSize',20)


figure(2);clf; hold on;
plot(V_ven_stack,P_ven_stack); 
plot(V_ven_avg,P_ven_avg,'k','LineWidth',3);
set(gca,'FontSize',20)



