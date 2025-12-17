function [time, ...
          sigma_LV, sigma_SEP, sigma_RV, ...
          Lsc_LV, Lsc_SEP, Lsc_RV, ...
          P_RV, P_PA, ...
          V_LV, V_RV, ...
          Cm_LV, Cm_SEP, Cm_RV, ...
          Q_a_valve, Q_m_valve] = load_outputs(animal_id,nx_or_hx_flag)

%{
    This function unpacks the "outputs" data structure. 
    Inputs: 
    animal_id     - animal number 
    nx_or_hx_flag - = 0 for nx and = 1 for hx
    Outputs: 
    Unpacked outputs
%}

% Access the data file, check if Nx or Hx 
    if nx_or_hx_flag == 0
        filename = sprintf('outputs_Nx%d.mat', animal_id);
        load(filename, 'outputs');
    elseif nx_or_hx_flag == 1
        filename = sprintf('outputs_Hx%d.mat', animal_id);
        load(filename, 'outputs');
    end

% Unpack the outputs 
    time = outputs.time;

    sigma_LV  = outputs.stresses.total.sigma_LV  * 7.5;
    sigma_SEP = outputs.stresses.total.sigma_SEP * 7.5;
    sigma_RV  = outputs.stresses.total.sigma_RV  * 7.5;

    Lsc_LV  = outputs.lengths.Lsc_LV;
    Lsc_SEP = outputs.lengths.Lsc_SEP;
    Lsc_RV  = outputs.lengths.Lsc_RV;

    P_RV = outputs.pressures.P_RV;
    P_PA = outputs.pressures.P_PA;

    V_LV = outputs.volumes.V_LV;
    V_RV = outputs.volumes.V_RV;

    Cm_LV  = outputs.curvatures.Cm_LV;
    Cm_SEP = outputs.curvatures.Cm_SEP;
    Cm_RV  = outputs.curvatures.Cm_RV;

    Q_a_valve = outputs.flows.Q_a_valve; 
    Q_m_valve = outputs.flows.Q_m_valve;

end