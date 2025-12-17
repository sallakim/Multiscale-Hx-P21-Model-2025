
addpath data_struct

data = makedatastructure_Nx_rat11;
HR = data.HR; 
T = 60/HR; 
dt = data.dt;
t = [0:dt:2*T]; 

para_fitted_Ca = [2	3	4	5	6	7	8	9	10;
0.0838	0.1306	0.1802	0.2557	0.3099	0.3613	0.408	0.4539	0.4602;
0.7513	0.8756	1.0274	1.4988	1.6107	1.6741	1.7902	2.1398	1.9832;
2.237	2.0486	1.948	1.852	1.6932	1.6773	1.5988	1.4952	1.4524;
0.1865	0.1815	0.1709	0.1693	0.161	0.1661	0.1425	0.1458	0.1222];
freq_all = para_fitted_Ca(1,:);
A_HR_pchip = pchip(freq_all,para_fitted_Ca(2,:));
A_HR = ppval(A_HR_pchip,HR/60);
B_HR_pchip = pchip(freq_all,para_fitted_Ca(3,:));
B_HR = ppval(B_HR_pchip,HR/60);
C_HR_pchip = pchip(freq_all,para_fitted_Ca(4,:));
C_HR = ppval(C_HR_pchip,HR/60);
Ca0_HR_pchip = pchip(freq_all,para_fitted_Ca(5,:));
Ca0_HR = ppval(Ca0_HR_pchip,HR/60);

a = A_HR; 
b = B_HR; 
c = C_HR; 
Ca0 = Ca0_HR;


stim_period = 1/(HR/60); 
phi = mod(t+0.0001,stim_period)/stim_period;
Ca_i = (a./phi).*exp(-b.*(log(phi)+c).^2) + Ca0;

%%
dt = 0.0001;
T = 1; 
t = [0:dt:2*T]; 
% Percentage of cardiac cycle 
k_TS = 0.04; % Beginning of cardiac cycle to maximal systole  
k_TR = 0.38; 

% Time to maximal systole 
TS_v = k_TS * T; 

% Time from maximal systole to relaxation 
TR_v = k_TR * T; 

eta = 2;
y = length(t);
for i = 1:length(t)
    tc = t(i);
    if tc >= 0 && tc < TS_v || tc >= T && tc < T + TS_v
    y(i) = 0.1610+ eta*0.5*(1 - cos(pi*(tc)/TS_v)); 
    elseif tc >= TS_v && tc < TR_v + TS_v || tc >= T+TS_v && tc < T+TR_v + TS_v
    y(i) = 0.1610+ eta*0.5*(1 + cos(pi*(tc - TS_v)/TR_v)); 
    elseif tc >= TR_v + TS_v && tc < T || tc >= T + TR_v + TS_v && tc < T+T
    y(i) = 0.1610; 
end 
end


hfig1 = figure(1);
clf
hold on
% plot (t,Ca_i,'--')
plot(t(1:10000),y(1:10000),'-')
% legend ('bahador','kim scale by x2')
print(hfig1,'-depsc2',strcat('Figures/','/Ca_Curve.eps'))
xticks([0,0.5,1])

