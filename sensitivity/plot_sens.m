function [] = plot_sens(y0,pars,sens,params,rat_id)

%% Plot 
sens = abs(sens);

% figure(5)
% clf
% plot(sens(1:329,:))
% title('V_{LV}')
% 
% figure(6)
% clf
% plot(sens(330:658,:))
% title('V_{RV}')
% 
% figure(7)
% clf
% plot(sens(659:987,:))
% title('P_{LV}')
% 
% figure(8)
% clf
% plot(sens(988:1316,:))
% title('R_{RV}')
% legend(params)

% Sens Scalar (ranked) 

% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');

% xlim([0 length(Isens)+1])
% xlabel('Parameters')
% Xlabel = params(Isens); 
% Xtick = 1:length(Isens);
% set(gca,'xtick',Xtick)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'XTickLabels',Xlabel)

% normalize sens to data
sens(1:329,:) = sens(1:329,:)/mean(y0(1:329));
sens(330:658,:) = sens(330:658,:)/mean(y0(330:658));
sens(659:987,:) = sens(659:987,:)/mean(y0(659:987));
sens(988:end,:) = sens(988:end,:)/mean(y0(988:end));

figure(100)
clf
hold on
sens_scalar_full =sqrt(sum(sens.^2)); 
[sens_scalar_full, sens_order_full] = sort(sens_scalar_full,'descend');
params_sorted_full = params(sens_order_full);
bar(sens_scalar_full./max(sens_scalar_full))
% plot(sens_scalar_full./max(sens_scalar_full),'k*')
xticks(1:length(pars))
xticklabels(params_sorted_full)
title(sprintf('Rat %d sens scalar',rat_id))
set(gca,'TickLabelInterpreter','latex')

figure(101)
clf
hold on
% plot(sens_scalar_full(1:10),'*')
bar(sens_scalar_full(1:10)./max(sens_scalar_full))
xticks(1:10)
xticklabels(params_sorted_full)
title(sprintf('Rat %d sens scalar: top 10',rat_id))
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')

figure(102)
clf
hold on

subplot(2,2,1)
colormap(cool)
sens_V_LV = sens(1:329,:);
sens_scalar =sqrt(sum(sens_V_LV.^2));
[sens_scalar, sens_order] = sort(sens_scalar,'descend');
params_sorted = params(sens_order);
% plot(sens_scalar,'*')
bar(sens_scalar(1:10),'FaceColor',[0.6,0.4,1.0],'EdgeColor',[0.6,0.4,1.0])
title('V_{LV}')
xticks(1:10)
xticklabels(params_sorted)
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')

top_params_vec = string(params_sorted(1:10)); 

subplot(2,2,2)
sens_V_LV = sens(330:658,:);
sens_scalar =sqrt(sum(sens_V_LV.^2));
[sens_scalar, sens_order] = sort(sens_scalar,'descend');
params_sorted = params(sens_order);
% plot(sens_scalar,'*')
bar(sens_scalar(1:10),'FaceColor',[0.6,0.4,1.0],'EdgeColor',[0.6,0.4,1.0])
title('V_{RV}')
xticks(1:10)
xticklabels(params_sorted)
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')

top_params_vec = [top_params_vec,string(params_sorted(1:10))];

subplot(2,2,3)
sens_V_LV = sens(659:987,:);
sens_scalar =sqrt(sum(sens_V_LV.^2));
[sens_scalar, sens_order] = sort(sens_scalar,'descend');
params_sorted = params(sens_order);
% plot(sens_scalar,'*')
bar(sens_scalar(1:10),'FaceColor',[0.6,0.4,1.0],'EdgeColor',[0.6,0.4,1.0])
title('P_{LV}')
xticks(1:10)
xticklabels(params_sorted)
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')

top_params_vec = [top_params_vec,string(params_sorted(1:10))];

subplot(2,2,4)
sens_V_LV = sens(988:end,:);
sens_scalar =sqrt(sum(sens_V_LV.^2));
[sens_scalar, sens_order] = sort(sens_scalar,'descend');
params_sorted = params(sens_order);
% plot(sens_scalar,'*')
bar(sens_scalar(1:10),'FaceColor',[0.6,0.4,1.0],'EdgeColor',[0.6,0.4,1.0])
title('P_{RV}')
xticks(1:10)
xticklabels(params_sorted)
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')

top_params_vec = [top_params_vec,string(params_sorted(1:10))];
top_params_vec = unique(top_params_vec,'stable');
top_params_vec = join(top_params_vec);

sgtitle(sprintf('Rat %d Most influential: \n%s',rat_id,top_params_vec),Interpreter='latex')


%% Condition number, corr matrix 
sens_top = sens(:,sens_order_full(1:6));

F = sens_top'*sens_top; 
F_cond = cond(F)

F_inv = inv(F); 

for i  = 1:1:6
    for j = 1:1:6
        corr_mat(i,j) = F_inv(i,j)/sqrt(F_inv(i,i)*F_inv(j,j));
    end
end
isupper = logical(triu(ones(size(corr_mat)),1)); 
corr_mat(isupper) = NaN; 
F_cond_rounded = round(F_cond);

% Plot
% clf
% figure(103)
% h = heatmap(abs(corr_mat),'MissingDataColor','w',Interpreter="latex");
% title(sprintf('Rat %s Correlation Matrix (absolute values), Condition #: %.0f',rat_label,(F_cond_rounded)));
% labels = params(sens_order_full(1:6));
% h.XDisplayLabels = labels;
% h.YDisplayLabels = labels; 

end