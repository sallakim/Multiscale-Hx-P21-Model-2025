function [] = plot_PV_loops(outputs,i_exp,label,nx_or_hx_flag,V_LV_stack,V_LV_avg,P_LV_stack,P_LV_avg,V_RV_stack,V_RV_avg,P_RV_stack,P_RV_avg)

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

minmax_LV = [min(min(V_LV_stack));max(max(V_LV_stack))];
minmax_RV = [min(min(V_RV_stack));max(max(V_RV_stack))];
A = [V_LV;V_RV;minmax_LV/1000;minmax_RV/1000];
minx = round(min(A)-0.001,3);
maxx = round(max(A)+0.001,3);

time = outputs.time;

if nx_or_hx_flag == 0
    hfig4a = figure(400); 
    sgtitle('RV: Optimized Model + Data')
    subplot(4,3,i_exp)
    hold on 
    % plot(V_RV_stack/1000,P_RV_stack,'Color',[.9 .9 .9],'LineWidth',1)
    plot(V_RV_avg/1000,P_RV_avg,'k','LineWidth',1)
    plot(V_RV,P_RV,'b','LineWidth',2)
    title(label)
    ylim([0 50])
    xlim([0,0.08])

    hfig4c = figure(402);
    sgtitle('LV: Optimized Model + Data')
    subplot(4,3,i_exp)
    hold on 
    % plot(V_LV_stack/1000,P_LV_stack,'Color',[.9 .9 .9],'LineWidth',1)
    plot(V_LV_avg/1000,P_LV_avg,'k','LineWidth',1)
    plot(V_LV,P_LV,'r','LineWidth',2)
    title(label)
    ylim([0 100])
    xlim([0,0.08])
end

if nx_or_hx_flag == 1
    hfig4b = figure(401); 
    sgtitle('RV: Optimized Model + Data')
    subplot(4,3,i_exp)
    hold on 
    % plot(V_RV_stack/1000,P_RV_stack,'Color',[.9 .9 .9],'LineWidth',1)
    plot(V_RV_avg/1000,P_RV_avg,'k','LineWidth',1)
    plot(V_RV,P_RV,'b','LineWidth',2)
    title(label)
    ylim([0 50])
    xlim([0,0.08])

    hfig4d = figure(403);
    sgtitle('LV: Optimized Model + Data')
    subplot(4,3,i_exp)
    hold on 
    % plot(V_LV_stack/1000,P_LV_stack,'Color',[.9 .9 .9],'LineWidth',1)
    plot(V_LV_avg/1000,P_LV_avg,'k','LineWidth',1)
    plot(V_LV,P_LV,'r','LineWidth',2)
    title(label)
    ylim([0 100])
    xlim([0,0.08])
end

end