% run optimization

clear all

% figure (1000)
% clf

% n_cores = 12;
% delete(gcp('nocreate'))
% p = parpool('local',n_cores);

% for animal = [11 12 51 52 54 55 56] 
for animal = [11] 
    clear data
    if animal==11
        % data = makedatastructure_Nx_rat11; % ---------------------------------------- change for each animal
        load ('Nx11_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==12
        % data = makedatastructure_Nx_rat12; % ---------------------------------------- change for each animal
        load ('Nx12_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==51
        % data = makedatastructure_Nx_rat51; % ---------------------------------------- change for each animal
        load ('Nx51_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==52
        % data = makedatastructure_Nx_rat52; % ---------------------------------------- change for each animal
        load ('Nx52_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==54
        % data = makedatastructure_Nx_rat54; % ---------------------------------------- change for each animal
        load ('Nx54_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==55
        % data = makedatastructure_Nx_rat55; % ---------------------------------------- change for each animal
        load ('Nx55_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    elseif animal==56
        % data = makedatastructure_Nx_rat56; % ---------------------------------------- change for each animal
        load ('Nx56_data.mat','V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean') % --- change for each animal
    end
    % table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
        table = readtable('P21_data_input.xlsx','PreserveVariableNames',true);
    data = make_datastructure_P21(animal,table);
    %% Change for each animal
    filename = strcat('opt_pars_Nx',num2str(animal),'.mat'); % ------------------------------------------- change for each animal


    %% Change for Nx vs Nx
    data.Hx_flag = 0; % ------------------------------------------ comment out for Nx
    [pars,UB,LB,data] = parameters_Nx(data); % ---------------------------------- change for each animal

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

    options = optimoptions('lsqnonlin','Display','iter','UseParallel',true,...
        'FiniteDifferenceStepSize',1e-3,'MaxFunctionEvaluations',2000);
    % xopt = lsqnonlin(fun,optx,optlow,opthi,options);

    % multistart
    optx_ms = ones(length(optx),21);
    optx_ms(:,1) = optx;
    size_optx_ms = size(optx_ms);
    i_length = size_optx_ms(2);
    j = 1;

    % initialize
    xopt_mat = ones(10,20);
    J_mat = ones(1,20);
    fun(optx);
    for i = 2:9%i_length
        opt_low_test = log( max(exp(optx).*0.5,exp(optlow)) );
        opt_hi_test = log( min(exp(optx).*1.5,exp(opthi)) );
        
        optx_ms(:,i) = unifrnd(opt_low_test,opt_hi_test);
    try
        [xopt,J] = lsqnonlin(fun,optx_ms(:,i),optlow,opthi,options);
        xopt_mat(:,j) = xopt;
        %     rout_mat(:,j) = rout;
        J_mat(j) = J;
        j = j+1;
    catch
        j=j+1;
    end
    end

    loc = find(J_mat == min(J_mat));
    xopt = xopt_mat(:,loc);

    %% Only want to change parameters in INDMAP

    pars_opt = pars;
    % pars_opt(INDMAP) = xopt;

    % save (filename)
end

elapsed_time = toc;
elapsed_time = elapsed_time/60

