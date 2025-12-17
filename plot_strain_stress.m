clear; clc; close all; 
addpath opt_pars\
addpath model\

nx_animal_ids = [11 12 51 52 54 55 56];
hx_animal_ids = [1 4 5 7 8 9 10 57 58 59 61 62]; 


%% Unpack outputs 

for i = 1:length(nx_animal_ids)
    animal_id = nx_animal_ids(i);
    filename = sprintf('outputs_Nx%d.mat', animal_id);
    load(filename)

    eps_LV = outputs.strains.eps_LV;
    eps_SEP = outputs.strains.eps_SEP; 
    eps_RV = outputs.strains.eps_RV; 

    Lsc_LV = outputs.lengths.Lsc_LV;
    Lsc_SEP = outputs.lengths.Lsc_SEP;
    Lsc_RV = outputs.lengths.Lsc_RV;
    
    sigma_pas_LV = outputs.stresses.passive.sigma_pas_LV;
    sigma_pas_SEP = outputs.stresses.passive.sigma_pas_SEP;
    sigma_pas_RV = outputs.stresses.passive.sigma_pas_RV;

    figure(1) 
    subplot(1,3,1)
    hold on
    shortening_LV = (Lsc_LV - max(Lsc_LV))./max(Lsc_LV);
    plot(shortening_LV,sigma_pas_LV)
    title('LV')
    xlabel('shortening')
    ylabel('stress')

    subplot(1,3,2)
    hold on 
    shortening_SEP = (Lsc_SEP - max(Lsc_SEP))./max(Lsc_SEP);
    plot(shortening_SEP,sigma_pas_SEP)
    title('SEP')
    xlabel('shortening')
    ylabel('stress')

    subplot(1,3,3)
    hold on 
    shortening_RV = (Lsc_RV - max(Lsc_RV))./max(Lsc_RV);
    plot(shortening_RV,sigma_pas_RV)
    title('RV')
    xlabel('shortening')
    ylabel('stress')

    sgtitle('Nx')

    dydx_max_nx(i) = max(diff(sigma_pas_RV)./diff(eps_RV));
    dydx_mean_nx(i) = mean(diff(sigma_pas_RV)./diff(eps_RV));
end


for i = 1:length(hx_animal_ids)
    animal_id = hx_animal_ids(i);
    filename = sprintf('outputs_Hx%d.mat', animal_id);
    load(filename)

    eps_LV = outputs.strains.eps_LV;
    eps_SEP = outputs.strains.eps_SEP; 
    eps_RV = outputs.strains.eps_RV; 

    Lsc_LV = outputs.lengths.Lsc_LV;
    Lsc_SEP = outputs.lengths.Lsc_SEP;
    Lsc_RV = outputs.lengths.Lsc_RV;
    
    sigma_pas_LV = outputs.stresses.passive.sigma_pas_LV;
    sigma_pas_SEP = outputs.stresses.passive.sigma_pas_SEP;
    sigma_pas_RV = outputs.stresses.passive.sigma_pas_RV;

    figure(2)
    subplot(1,3,1)
    hold on 
    shortening_LV = (Lsc_LV - max(Lsc_LV))./max(Lsc_LV);
    plot(shortening_LV,sigma_pas_LV)
    title('LV')
    xlabel('shortening')
    ylabel('stress')

    subplot(1,3,2)
    hold on 
    shortening_SEP = (Lsc_SEP - max(Lsc_SEP))./max(Lsc_SEP);
    plot(shortening_SEP,sigma_pas_SEP)
    title('SEP')
    xlabel('shortening')
    ylabel('stress')

    subplot(1,3,3)
    hold on 
    shortening_RV = (Lsc_RV - max(Lsc_RV))./max(Lsc_RV);
    plot(shortening_RV,sigma_pas_RV)
    title('RV')
    xlabel('shortening')
    ylabel('stress')

    sgtitle('Hx')

    dydx_max_hx(i) = max(diff(sigma_pas_RV)./diff(eps_RV));
    dydx_mean_hx(i) = mean(diff(sigma_pas_RV)./diff(eps_RV));

end
