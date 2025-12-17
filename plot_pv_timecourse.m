function [] = plot_pv_timecourse(outputs,i_exp,label,nx_or_hx_flag,V_LV_stack,V_LV_avg,P_LV_stack,P_LV_avg,V_RV_stack,V_RV_avg,P_RV_stack,P_RV_avg)

%{
This function plots the model PV loops with the data. 
inputs: 
"outputs" - model solutions 
"i_exp"   - experiment index 
"animal_id" - current animal number 
"nx_ids"    - vector of the Nx animal ids 
"V_i_stack" - volume data for i = LV, RV 
"P_i_stack" - pressure data 
"V_i_avg"   - average volume data 
"P_i_avg"   - average pressure data 
outputs: 
figure with the solved model PV loop and PV loop data 
%}

%% Unpack outputs 
V_LV = outputs.volumes.V_LV; 
V_RV = outputs.volumes.V_RV; 

P_LV = outputs.pressures.P_LV; 
P_RV = outputs.pressures.P_RV; 

P_PA = outputs.pressures.P_PA; 
P_SA = outputs.pressures.P_SA; 

minmax_LV = [min(min(V_LV_stack));max(max(V_LV_stack))];
minmax_RV = [min(min(V_RV_stack));max(max(V_RV_stack))];
A = [V_LV;V_RV;minmax_LV/1000;minmax_RV/1000];
minx = round(min(A)-0.001,3);
maxx = round(max(A)+0.001,3);

time = outputs.time;
%
loc_half_time = round(length(time)/2); 
time_50 = linspace(0,time(loc_half_time),50);

grey = [.8 .8 .8]; 

if nx_or_hx_flag == 0
    hfig4a = figure(405); 
    sgtitle('RV: Pressure')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,P_RV_stack,'Color',grey,'LineWidth',1)
    plot(time_50,P_RV_avg,'k','LineWidth',1)
    plot(time(1:loc_half_time),P_RV(1:loc_half_time),'b','LineWidth',2)
    title('Nx',label)

    hfig4b = figure(406); 
    sgtitle('RV: Volume')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,V_RV_stack/1000,'Color',grey,'LineWidth',1)
    plot(time_50,V_RV_avg/1000,'k','LineWidth',1)
    plot(time(1:loc_half_time),V_RV(1:loc_half_time),'b','LineWidth',2)
    title('Nx',label)

    hfig4c = figure(407); 
    sgtitle('LV: Pressure')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,P_LV_stack,'Color',grey,'LineWidth',1)
    plot(time_50,P_LV_avg,'k','LineWidth',1)
    plot(time(1:loc_half_time),P_LV(1:loc_half_time),'r','LineWidth',2)
    title('Nx',label)

    hfig4d = figure(408); 
    sgtitle('LV: Volume')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,V_LV_stack/1000,'Color',grey,'LineWidth',1)
    plot(time_50,V_LV_avg/1000,'k','LineWidth',1)
    plot(time(1:loc_half_time),V_LV(1:loc_half_time),'r','LineWidth',2)
    title('Nx',label)

    hfig4a = figure(413); 
    sgtitle('PA and Ao: Pressure')
    subplot(3,4,i_exp)
    hold on 
    plot(time(1:loc_half_time),P_PA(1:loc_half_time),'b:','LineWidth',2)
    plot(time(1:loc_half_time),P_SA(1:loc_half_time),'r','LineWidth',2)
    title('Nx',label)

end

if nx_or_hx_flag == 1
    hfig4e = figure(409); 
    sgtitle('RV: Pressure')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,P_RV_stack,'Color',grey,'LineWidth',1)
    plot(time_50,P_RV_avg,'k','LineWidth',1)
    plot(time(1:loc_half_time),P_RV(1:loc_half_time),'b','LineWidth',2)
    title('Hx',label)

    hfig4f = figure(410); 
    sgtitle('RV: Volume')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,V_RV_stack/1000,'Color',grey,'LineWidth',1)
    plot(time_50,V_RV_avg/1000,'k','LineWidth',1)
    plot(time(1:loc_half_time),V_RV(1:loc_half_time),'b','LineWidth',2)
    title('Hx',label)

    hfig4g = figure(411); 
    sgtitle('LV: Pressure')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,P_LV_stack,'Color',grey,'LineWidth',1)
    plot(time_50,P_LV_avg,'k','LineWidth',1)
    plot(time(1:loc_half_time),P_LV(1:loc_half_time),'r','LineWidth',2)
    title('Hx',label)

    hfig4h = figure(412); 
    sgtitle('LV: Volume')
    subplot(3,4,i_exp)
    hold on 
    plot(time_50,V_LV_stack/1000,'Color',grey,'LineWidth',1)
    plot(time_50,V_LV_avg/1000,'k','LineWidth',1)
    plot(time(1:loc_half_time),V_LV(1:loc_half_time),'r','LineWidth',2)
    title('Hx',label)

    hfig4a = figure(415); 
    sgtitle('PA and Ao: Pressure')
    subplot(3,4,i_exp)
    hold on 
    plot(time(1:loc_half_time),P_PA(1:loc_half_time),'b:','LineWidth',2)
    plot(time(1:loc_half_time),P_SA(1:loc_half_time),'r','LineWidth',2)
    title('Hx',label)

end

end