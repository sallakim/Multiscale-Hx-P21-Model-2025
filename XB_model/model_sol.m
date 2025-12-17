function [outputs] = model_sol(data,pars,init)

%% Initialization  

tspan = data.tspan;  
dt    = data.dt; 

T = data.T; 

ODE_TOL = data.gpars.ODE_TOL;


%% Set mass matrix M for DAE 
M = speye(length(init));
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 

%% Solve model 

% Solve for 20 beats 
opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);
sol  = ode15s(@model_xb,[tspan(1) tspan(end)],init,opts,pars,data);
sols = deval(sol,tspan);  

t = sol.x(1):dt:sol.x(end); 
beats = mod(t,T);
x = find(round(beats,3) == 0);
y = find(t(x) <= t(end)); 
init  = sols(:,x(y(end-1))); 

% Solve 1 more beat
time = [0:dt:T]; 
sol  = ode15s(@model_xb,[time(1) time(end)],init,opts,pars,data);
sols = deval(sol,time);  % sols length (# initial conditions, # of timesteps) 
sols = sols'; 

%% Calculate other time-varying model quantities (pressures, flows, etc.) 
for i = 1:length(time) 
    [~,o(:,i)] = model_xb(time(i),sols(i,:),pars,data);
end 

%% Create output structure  
outputs.stressess.sigma_act_LV = o(1,:); 
outputs.stressess.sigma_act_RV = o(3,:); 
outputs.stressess.sigma_RV = o(6,:); 
outputs.stressess.sigma_pas_RV = o(9,:); 
outputs.time = time;

end