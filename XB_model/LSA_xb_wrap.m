% LSA Model Wrap 

function y = LSA_xb_wrap(pars,data) 

outputs = model_sol(pars,data);

sigma_act_RV = outputs.stressess.sigma_act_RV;

y = [sigma_act_RV];

end