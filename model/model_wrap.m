%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function [rout,sol] = model_wrap(pars,data)

ALLPARS = data.gpars.ALLPARS;
INDMAP  = data.gpars.INDMAP;

tpars = ALLPARS;
tpars(INDMAP') = pars;

[sol,rout,J] = model_sol(tpars,data);

% figure (1)
% plot(rout)

% only want to load in the most significant parameters, reassign