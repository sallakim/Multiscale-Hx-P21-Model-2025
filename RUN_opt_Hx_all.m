% run optimization

clear all

% figure (1000)
clf

% n_cores = 12;
% delete(gcp('nocreate'))
% p = parpool('local',n_cores);

for animal = [1 4 5 7 8 9 10 57 58 59 61 62]
    clear data
    if animal==1
        % data = makedatastructure_Hx_rat1; % ---------------------------------------- change for each animal
        load ('new_data/Hx1_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==4
        % data = makedatastructure_Hx_rat4; % ---------------------------------------- change for each animal
        load ('new_data/Hx4_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==5
        % data = makedatastructure_Hx_rat5; % ---------------------------------------- change for each animal
        load ('new_data/Hx5_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==7
        % data = makedatastructure_Hx_rat7; % ---------------------------------------- change for each animal
        load ('new_data/Hx7_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==8
        % data = makedatastructure_Hx_rat8; % ---------------------------------------- change for each animal
        load ('new_data/Hx8_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==9
        % data = makedatastructure_Hx_rat9; % ---------------------------------------- change for each animal
        load ('new_data/Hx9_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==10
        % data = makedatastructure_Hx_rat10; % ---------------------------------------- change for each animal
        load ('new_data/Hx10_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==57
        % data = makedatastructure_Hx_rat57; % ---------------------------------------- change for each animal
        load ('new_data/Hx57_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==58
        % data = makedatastructure_Hx_rat58; % ---------------------------------------- change for each animal
        load ('new_data/Hx58_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==59
        % data = makedatastructure_Hx_rat59; % ---------------------------------------- change for each animal
        load ('new_data/Hx59_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==61
        % data = makedatastructure_Hx_rat61; % ---------------------------------------- change for each animal
        load ('new_data/Hx61_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==62
        % data = makedatastructure_Hx_rat62; % ---------------------------------------- change for each animal
        load ('new_data/Hx62_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    end
     table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
    data = make_datastructure_P21(animal,table);
    %% Change for each animal
    filename = strcat('opt_pars_Hx',num2str(animal),'_NEW_largebounds_SALLAUPDATE_2025.mat'); % ------------------------------------------- change for each animal


    %% Change for Nx vs Hx
    data.hx_flag = 1; % ------------------------------------------ comment out for Nx
    [pars,UB,LB,data] = parameters_Hx(data); % ---------------------------------- change for each animal

    %%
    INDMAP = [5, 6, 11, 12, 13, 14, 15, 16, 17, 18] % top 10


    ALLPARS  = pars;
    ODE_TOL  = 1e-6;
    DIFF_INC = sqrt(ODE_TOL);

    gpars.INDMAP   = INDMAP;
    gpars.ALLPARS  = ALLPARS;
    gpars.ODE_TOL  = ODE_TOL;
    gpars.DIFF_INC = DIFF_INC;

    data.gpars = gpars;

    data.V_LV_avg = V_LV_mean;
    data.V_RV_avg = V_RV_mean;
    data.P_LV_avg = P_LV_mean;
    data.P_RV_avg = P_RV_mean;

    data.MgATP_cytoplasm = 8;
    data.MgADP_cytoplasm = 0.05;
    data.Pi_cyto         = 1.3;


    %% Optimization

    optx   = pars(INDMAP);
    opthi  = UB(INDMAP);
    optlow = LB(INDMAP);
    opthi(end-1:end) = [log(0.1) log(0.5)];
    optlow(end-1:end) = [log(0.001) log(0.3)];
    


    fun = @(x) model_wrap(x,data);

    % options = optimoptions('fmincon','Display','iter','UseParallel',false);
    %xopt = fmincon(fun,optx,[],[],[],[],optlow,opthi,[],options);
    % options = optimoptions('lsqnonlin','Display','iter','UseParallel',true);

    options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','iter','UseParallel',true,'FiniteDifferenceStepSize',1e-3,'MaxFunctionEvaluations',2000);
    % xopt = lsqnonlin(fun,optx,optlow,opthi,options);

    % multistart
    optx_ms = ones(length(optx),21);
    optx_ms(:,1) = optx;
    size_optx_ms = size(optx_ms);
    i_length = size_optx_ms(2);
    j = 1;

    % initialize
    xopt_mat = ones(10,20);
    J_mat = ones(1,20).*1e6;

    for i = 2:9%i_length

        opt_low_test = log( max(exp(optx).*0.5,exp(optlow)) );
        opt_hi_test = log( min(exp(optx).*1.5,exp(opthi)) );
        
        optx_ms(:,i) = unifrnd(opt_low_test,opt_hi_test);
        % optx_ms(:,i) = unifrnd(min(optx.*0.5,optx.*1.5),max(optx.*0.5,optx.*1.5));
        % optx_ms(end-1:end,i) = log(unifrnd([0.01 0.01],[0.4,0.5]));
    try
        [xopt,J] = lsqnonlin(fun,optx_ms(:,i),optlow,opthi,options);
        xopt_mat(:,j) = xopt;
        %     rout_mat(:,j) = rout;
        J_mat(j) = J;
        j = j+1;

        [rout,sol] = model_wrap(xopt,data);
        % subplot(1,2,2);hold on; plot(sol.P_LV_model);




        figure(i); clf
        subplot(3,2,1); hold on;
        plot(P_LV_mean,'k'); plot(sol.P_LV_model,'r');
        subplot(3,2,2); hold on;
        plot(P_RV_mean,'k'); plot(sol.P_RV_model,'b');
        subplot(3,2,3); hold on;
        plot(V_LV_mean,'k'); plot(sol.V_LV_model.*1000,'r');
        subplot(3,2,4); hold on;
        plot(V_RV_mean,'k'); plot(sol.V_RV_model.*1000,'b');
        subplot(3,2,5); hold on;
        plot(V_LV_mean,P_LV_mean,'k'); plot(sol.V_LV_model.*1000,sol.P_LV_model,'r');
        subplot(3,2,6); hold on;
        plot(V_RV_mean,P_RV_mean,'k'); plot(sol.V_RV_model.*1000,sol.P_RV_model,'b');
    catch
        j=j+1;
    end
    end

    loc = find(J_mat == min(J_mat));
    xopt = xopt_mat(:,loc);

    %% Only want to change parameters in INDMAP

    pars_opt = pars;
    % pars_opt(INDMAP) = xopt;

    save (filename)
end

elapsed_time = toc;
elapsed_time = elapsed_time/60

