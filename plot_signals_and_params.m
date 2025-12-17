load opt_pars_Nx12.mat ALLPARS xopt_mat J_mat data P_LV_mean P_RV_mean V_LV_mean V_RV_mean
param_names =  {'C_{SA}', 'C_{SV}', 'C_{PA}', 'C_{PV}', ...
    'R_{SA}', 'R_{PA}', ...
    'R_m', 'R_a', 'R_t', 'R_p', ...
    'Amref_{LV}', 'Amref_{SEP}', 'Amref_{RV}', ...
    'Vw_{LV}', 'Vw_{SEP}', 'Vw_{RV}', ...
    'k_{TS}', 'k_{TR}'};
INDMAP = [5, 6, 11, 12, 13, 14, 15, 16, 17, 18];

addpath model\
n_starts = 20;
close all
[who,where] = sort(J_mat,'ascend');
counter = 1;
for i=where(1:20)%1:20
    try
        % figure(99);
        % [rout,sol] = model_wrap(optx_ms(:,i),data);
        % subplot(1,2,1);hold on; plot(sol.P_LV_model);
        [rout,sol] = model_wrap(xopt_mat(:,i),data);
        % subplot(1,2,2);hold on; plot(sol.P_LV_model);




        figure(1);
        subplot(2,2,1); hold on;
        plot(sol.pressures.P_LV,'LineWidth',1); ylabel('LV Pressure');
        set(gca,'FontSize',15); grid on; xticklabels({'0','T'})
        subplot(2,2,2); hold on;
        plot(sol.pressures.P_RV,'LineWidth',1); ylabel('RV Pressure');
        set(gca,'FontSize',15); grid on; xticklabels({'0','T'})
        subplot(2,2,3); hold on;
        plot(sol.volumes.V_LV.*1000,'LineWidth',1); ylabel('LV Volume');
        set(gca,'FontSize',15); grid on; xticklabels({'0','T'})
        subplot(2,2,4); hold on;
        plot(sol.volumes.V_RV.*1000,'LineWidth',1); ylabel('RV Volume');
        set(gca,'FontSize',15); grid on; xticklabels({'0','T'})
        


        figure(2); hold on;
        plot(xopt_mat(:,i),'o','MarkerSize',8,'LineWidth',2);
        grid on;
        set(gca,'FontSize',15); grid on;
        ylabel('Log Parameter Estimate');
        xlabel('Parameter');
        xticks(1:10);
        xticklabels(param_names(INDMAP))

        if counter>1
            figure(3);
            subplot(2,2,1); hold on;
            plot(PLV_best-sol.pressures.P_LV,'LineWidth',1); ylabel('LV Pressure');
            set(gca,'FontSize',15); grid on;
            subplot(2,2,2); hold on;
            plot(PRV_best-sol.pressures.P_RV,'LineWidth',1); ylabel('RV Pressure');
            set(gca,'FontSize',15); grid on;
            subplot(2,2,3); hold on;
            plot(VLV_best-sol.volumes.V_LV.*1000,'LineWidth',1); ylabel('LV Volume');
            set(gca,'FontSize',15); grid on;
            subplot(2,2,4); hold on;
            plot(VRV_best-sol.volumes.V_RV.*1000,'LineWidth',1); ylabel('RV Volume');
            set(gca,'FontSize',15); grid on;
        else
            PLV_best = sol.pressures.P_LV;
            PRV_best = sol.pressures.P_RV;
            VLV_best = sol.volumes.V_LV.*1000;
            VRV_best = sol.volumes.V_RV.*1000;

        end

    catch
    end

    counter = counter+1;

end
