%% Analysis of metabolic features of quiescence deepening with fluxes 
% generated via CFR / linear iMAT

%% Initializing COBRA

clear;

initCobraToolbox;
changeCobraSolver('gurobi');


%% Analysis of flux-based model

% Loads model and fluxes
load QDS_model_3.mat


%% Correlation analysis
% which reactions correlate best with quiescence depth

corr_r = zeros([length(model.rxns) 1]);
corr_p = zeros([length(model.rxns) 1]);
corr_q = zeros([length(model.rxns) 1]);

for i = 1:length(model.rxns)

    [corr_r(i), corr_p(i)] = corr(flux_q(:, i), qds_flux);

end

[h, ~, ~, corr_q] = fdr_bh(corr_p);


%% Analysis of subsystems correlation coefficients

% replace NaN correlations with 0
corr_r(isnan(corr_r)) = 0;

% find total subsystems in all reactions
subsystems_c_m = model.subSystems;
subsystems_m = string();

for i = 1:length(subsystems_c_m)

    subsystems_1_m = (subsystems_c_m{i});
    subsystems_m(i) = string(subsystems_1_m(1, 1));

end

% aggregate reaction information: name, subsystem, correlation with
% quiescence deepening (pearson correlation coefficient)
lr_model_t = table(model.rxnNames, subsystems_m', corr_r, abs(corr_r), 'VariableNames', {'Reaction_names', 'Subsystems','Pearson_R', 'Pearson_R_Abs'});

% finds average pearson correlation coefficient for each subsystem
total_subsystems_unique = unique(subsystems_m);

coeff_avg = zeros([length(total_subsystems_unique) 1]); 

coeff_CI = zeros([length(total_subsystems_unique) 2]); 

t_test_pval = zeros([length(total_subsystems_unique) 1]); 

count = zeros([length(total_subsystems_unique) 1]); 

for i = 1:length(total_subsystems_unique)

    [coeff_avg(i), ~, coeff_CI(i, :)] = normfit(corr_r(strcmp(lr_model_t.Subsystems, total_subsystems_unique(i))));
    
    [~, t_test_pval(i)] = ttest(corr_r(strcmp(lr_model_t.Subsystems, total_subsystems_unique(i))));
    
    count(i) = length(corr_r(strcmp(lr_model_t.Subsystems, total_subsystems_unique(i))));

end


% benjamini-hochberg FDR p-value correction
[~, ~, ~, qval] = fdr_bh(t_test_pval);

% creates table for subsystems
subsystem_lr_model_t = table(total_subsystems_unique', coeff_avg, coeff_CI(:, 1), coeff_CI(:, 2), qval, count,...
    'VariableNames', {'Subsystems', 'Mean_coeff_5th_root', 'CI_95_low', 'CI_95_high','q_value', 'count'});

subsystem_lr_model_t = sortrows(subsystem_lr_model_t, 2, 'descend');

% retains only significant subsystems (different from zero)
subsystem_lr_model_t = subsystem_lr_model_t(subsystem_lr_model_t.q_value < 5E-2, :)

