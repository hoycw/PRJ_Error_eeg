%% Model Predictions
cond_lab = {'Likely Win', 'Unlikely Neu', 'Unlikely Loss', 'Unlikely Win', 'Unlikely Neu', 'Likely Loss'};
% cond_lab = {'Easy Win', 'Easy Neu', 'Easy Loss', 'Hard Win', 'Hard Neu', 'Hard Loss'};

% Categorical Reward
exp_value   = [1 1 1 -1 -1 -1];
valence_cat   = [1 nan 0 1 nan 0];
valence       = [1 0 -1 1 0 -1];
magnitude   = [1 0 1 1 0 1];
interaction = exp_value .* valence;

% Continuous Reward Prediction Error
rp      = [0.6 0.6 0.6 -0.6 -0.6 -0.6];
rpe     = valence - rp;
rpe_mag = abs(rpe);

% Categorical Likelihood
lik_cat   = [1 0 0 0 0 1];
% Continuous Likelihood
lik_cont  = [0.75 0.1 0.15 0.15 0.1 0.75];
rare_cont = 1-lik_cont;

%% Plot Models
fig_name  = 'model_prediction_comparison';
fig_ftype = 'png';
figure('Name',fig_name,'units','normalized','OuterPosition',[0 0 1 1]);

valence_color = [209 151 105]./255;
valence_style = '-';
valence_mark  = 's';

mag_color = [118 160 156]./255;
mag_style = ':';
mag_mark  = 'o';

lik_color = [152 78 163]./255;
lik_style = '--';
lik_mark  = '*';

alt_color = 'k';    % categorical salience/likelihood/frequency
odd_color = 'r';    % meaningless interaction
alt_style = '-.';
alt_mark  = 's';

cond_idx = 1:6;
easy_idx = 1:3;
hard_idx = 4:6;

cond_idx_nn = [1 3 4 6];
easy_idx_nn = [1 3];
hard_idx_nn = [4 6];

offset     = 0.08;
tick_angle = -20;
width      = 2;

%% ======================= Binary Outcomes =======================
%--------------------- Binary Outcomes: Categorical Coding ---------------------
subplot(2,2,1); hold on;
title('Binary Outcomes: Categorical Coding');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx_nn);
set(gca,'XTickLabels',cond_lab(cond_idx_nn));
xtickangle(tick_angle);
xlim([0 7]);
ylim([-0.3 1.1]);
yticks([0 1]);

tmp1 = plot(easy_idx_nn, valence_cat(easy_idx_nn), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);
plot(hard_idx_nn, valence_cat(hard_idx_nn), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);

tmp2 = plot(easy_idx_nn+offset, magnitude(easy_idx_nn), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);
plot(hard_idx_nn+offset, magnitude(hard_idx_nn), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);

tmp3 = plot(easy_idx_nn+2*offset, lik_cat(easy_idx_nn), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
plot(hard_idx_nn+2*offset, lik_cat(hard_idx_nn), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);

% tmp4 = plot(easy_idx_nn+3*offset, interaction(easy_idx_nn), 'Color', alt_color, 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
% plot(hard_idx_nn+3*offset, interaction(hard_idx_nn), 'Color', alt_color, 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);

legend([tmp1 tmp2 tmp3], {'Reward Valence','Reward Magnitude','Reward Likelihood'}, 'Location','southwest');
% legend([tmp1 tmp2 tmp3 tmp4], {'Valence','Magnitude','Likelihood','Interaction'}, 'Location','southwest');
[~] = fn_min_white_space(gca);

%--------------------- Binary Outcomes: Prediction Errors ---------------------
subplot(2,2,2); hold on;
title('Binary Outcomes: Prediction Errors');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx_nn);
set(gca,'XTickLabels',cond_lab(cond_idx_nn));
xtickangle(tick_angle);
xlim([0 7]);
ylim([-2.1 2.1]);
yticks([-2:1:2]);

tmp1 = plot(easy_idx_nn, rpe(easy_idx_nn), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);
plot(hard_idx_nn, rpe(hard_idx_nn), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);

tmp2 = plot(easy_idx_nn+offset, rpe_mag(easy_idx_nn), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);
plot(hard_idx_nn+offset, rpe_mag(hard_idx_nn), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);

tmp3 = plot(easy_idx_nn+2*offset, lik_cont(easy_idx_nn), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
plot(hard_idx_nn+2*offset, lik_cont(hard_idx_nn), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
legend([tmp1 tmp2 tmp3], {'RPE','RPE Salience','Outcome Likelihood'}, 'Location','southeast');
[~] = fn_min_white_space(gca);

%% ======================= All Outcomes =======================
%--------------------- All Outcomes: Directional Coding ---------------------
subplot(2,2,3); hold on;
title('All Outcomes: Directional Coding');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx);
set(gca,'XTickLabels',cond_lab);
xtickangle(tick_angle);
xlim([0 7]);
ylim([-1.6 1.1]);
yticks([-1 0 1]);

tmp1 = plot(easy_idx, valence(easy_idx), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);
plot(hard_idx, valence(hard_idx), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);

tmp2 = plot(easy_idx+offset, magnitude(easy_idx), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);
plot(hard_idx+offset, magnitude(hard_idx), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);

tmp3 = plot(easy_idx+2*offset, lik_cat(easy_idx), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
plot(hard_idx+2*offset, lik_cat(hard_idx), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);

% tmp4 = plot(easy_idx+3*offset, interaction(easy_idx), 'Color', 'r', 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
% plot(hard_idx+3*offset, interaction(hard_idx), 'Color', 'r', 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);

legend([tmp1 tmp2 tmp3], {'Reward Valence','Reward Magnitude','Reward Likelihood'}, 'Location','southwest');
% legend([tmp1 tmp2 tmp3 tmp4], {'Valence','Magnitude','Likelihood', 'Interaction'}, 'Location','southwest');
[~] = fn_min_white_space(gca);

%--------------------- All Outcomes: Prediction Errors ---------------------
subplot(2,2,4); hold on;
title('All Outcomes: Prediction Errors');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx);
set(gca,'XTickLabels',cond_lab);
xtickangle(tick_angle);
xlim([0 7]);
ylim([-2.1 2.1]);
yticks([-2:1:2]);

tmp1 = plot(easy_idx, rpe(easy_idx), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);
plot(hard_idx, rpe(hard_idx), 'Color', valence_color, 'LineStyle', valence_style, 'LineWidth', width, 'Marker', valence_mark);

tmp2 = plot(easy_idx+offset, rpe_mag(easy_idx), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);
plot(hard_idx+offset, rpe_mag(hard_idx), 'Color', mag_color, 'LineStyle', mag_style, 'LineWidth', width, 'Marker', mag_mark);

tmp3 = plot(easy_idx+2*offset, lik_cont(easy_idx), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
plot(hard_idx+2*offset, lik_cont(hard_idx), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
legend([tmp1 tmp2 tmp3], {'RPE','RPE Salience','Outcome Likelihood'}, 'Location','southeast');
[~] = fn_min_white_space(gca);

%% Save figure
fig_dir = [root_dir 'PRJ_Error_eeg/results/model_predictions/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

fig_fname = [fig_dir fig_name '.' fig_ftype];
% % Commented out because screen ratios aren't right automatically
fprintf('Saving %s\n',fig_fname);
saveas(gcf,fig_fname);
