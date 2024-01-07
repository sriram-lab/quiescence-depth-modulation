 %% Initial analysis of metabolic changes associated with quiescence deepening

% Analysis of metabolic fluxes generated through INIT + CHRR

%% Initializing Cobra + Model

clear;

% initializes cobra toolbox
initCobraToolbox;
changeCobraSolver('gurobi');

%% Loads model and fluxes

load Rat-GEM.mat
model = ratGEM;
model.b = model.b(:, 1);
model = generateRules(model);
model.genes = upper(model.genes);
model.c(model.c == 1) = 0;

load init_chrr_alpha_1.1_beta_0.6_all_qds.MAT


%% Correlation analysis
% which reactions correlate best with quiescence depth

flux_init_chrr = cell2mat(flux_avg_array);

qds = [2 3 4 6 8 10 12 14 16];
corr_r = zeros([length(model.rxns) 1]);
corr_p = zeros([length(model.rxns) 1]);
corr_q = zeros([length(model.rxns) 1]);

for i = 1:length(model.rxns)

    [corr_r(i), corr_p(i)] = corr(flux_init_chrr(i, :)', qds');

end

[~, ~, ~, corr_q] = fdr_bh(corr_p);

%% top correlating reactions

reactions_corr_table = table(model.rxns, corr_r, corr_p, corr_q, 'VariableNames', {'rxn', 'r', 'p', 'q'})

reactions_corr_table = sortrows(reactions_corr_table, 2, 'descend')

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


%% Analysis of biological significance: subsystem analysis

%% Up in quiescence deepening
% enrichment analysis:

% find subsystems in subset of reactions and converts to 1x1 string

% up in quiescence deepening
subsystems_c = (model.subSystems(corr_r > 0.0 & corr_p < 0.15))

% all
% subsystems_c = (model.subSystems(corr_q < 0.05))

subsystems_subset = string();

for i = 1:length(subsystems_c)

    subsystems_1 = (subsystems_c{i});
    subsystems_subset(i) = string(subsystems_1(1, 1));

end

% find total subsystems in all reactions and converts to 1x1 string
subsystems_c_m = model.subSystems;
subsystems_m = string();

for i = 1:length(subsystems_c_m)
i;
    subsystems_1_m = (subsystems_c_m{i});
    subsystems_m(i) = string(subsystems_1_m(1, 1));

end

subsystem_types = unique(subsystems_subset);
fold_enrichment = zeros([length(subsystem_types) 1]);
pval = zeros([length(subsystem_types) 1]);

%
for i = 1:length(subsystem_types)

    subsystem_types(i);

    num_sample = sum(strcmp(subsystems_subset, subsystem_types(i)));

    num_overall = sum(strcmp(subsystems_m, subsystem_types(i)));

    fold_enrichment(i) = (num_sample ./ length(subsystems_subset)) ./ (num_overall ./ length(model.subSystems));

    % enrichment p-value calculation
    pval(i) = 1 - hygecdf(num_sample, length(subsystems_m), num_overall, length(subsystems_subset));

end

% benjamini-hochberg FDR p-value correction
[h, p_crit, ~, qval] = fdr_bh(pval, 0.05);


pathway_analysis_table = table(subsystem_types', fold_enrichment, qval, 'VariableNames', {'Subsystem', 'Fold_enrichment', 'Qval'});

% sort by fold enrichment
pathway_analysis_table = sortrows(pathway_analysis_table, 2, 'descend')

% retain only q-value < 0.01
pathway_analysis_table = pathway_analysis_table(pathway_analysis_table.Qval < 0.01, :)


%% Down in quiescence deepening
% enrichment analysis:

% find subsystems in subset of reactions and converts to 1x1 string

% down in quiescence deepening
subsystems_c = (model.subSystems(corr_r < 0.0 & corr_p < 0.15))

% all
% subsystems_c = (model.subSystems(corr_q < 0.05))

subsystems_subset = string();

for i = 1:length(subsystems_c)

    subsystems_1 = (subsystems_c{i});
    subsystems_subset(i) = string(subsystems_1(1, 1));

end

% find total subsystems in all reactions and converts to 1x1 string
subsystems_c_m = model.subSystems;
subsystems_m = string();

for i = 1:length(subsystems_c_m)
i;
    subsystems_1_m = (subsystems_c_m{i});
    subsystems_m(i) = string(subsystems_1_m(1, 1));

end

subsystem_types = unique(subsystems_subset);
fold_enrichment = zeros([length(subsystem_types) 1]);
pval = zeros([length(subsystem_types) 1]);


for i = 1:length(subsystem_types)

    subsystem_types(i);

    num_sample = sum(strcmp(subsystems_subset, subsystem_types(i)));

    num_overall = sum(strcmp(subsystems_m, subsystem_types(i)));

    fold_enrichment(i) = (num_sample ./ length(subsystems_subset)) ./ (num_overall ./ length(model.subSystems));

    % enrichment p-value calculation
    pval(i) = 1 - hygecdf(num_sample, length(subsystems_m), num_overall, length(subsystems_subset));

end

% benjamini-hochberg FDR p-value correction
[h, p_crit, ~, qval] = fdr_bh(pval, 0.05);

pathway_analysis_table = table(subsystem_types', fold_enrichment, qval, 'VariableNames', {'Subsystem', 'Fold_enrichment', 'Qval'});

% sort by fold enrichment
pathway_analysis_table = sortrows(pathway_analysis_table, 2, 'descend')

% retain only q-value < 0.01
pathway_analysis_table = pathway_analysis_table(pathway_analysis_table.Qval < 0.01, :)


%% time vs flux for positively and negatively correlating reactions

flux_init_chrr_z = normalize(flux_init_chrr, 2)

%% up 
up_rxns_quies_deepening = corr_r > 0 & corr_p < 0.15;

sum(up_rxns_quies_deepening)

plot(qds, flux_init_chrr_z(up_rxns_quies_deepening, :),...
    'LineWidth', 0.4,...
    'Color',[0.0824 0.3765 0.7412 0.3])
xlabel("Quiescence depth (days)")
ylabel("Average flux (mmol / (gDW hr), z-score)")


%% down
dw_rxns_quies_deepening = corr_r < 0 & corr_p < 0.15;

sum(dw_rxns_quies_deepening)

plot(qds, flux_init_chrr_z(dw_rxns_quies_deepening, :),...
    'LineWidth', 0.4,...
    'Color',[0.0824 0.3765 0.7412 0.4])
xlabel("Quiescence depth (days)")
ylabel("Average flux (mmol / (gDW hr), z-score)")


%% clustering of transcriptomics used in model reduction 

% Loads quiescence deepening transcriptomics
fujimaki_19 = readtable('fujimaki_19_GSE124109.txt');

% isolates genes in model
fujimaki_19.Gene = upper(fujimaki_19.Gene)
fujimaki_19 = fujimaki_19(ismember(fujimaki_19.Gene, model.genes), :)

% create matrix + omit controls
fujimaki_19_mat = table2array(fujimaki_19(:, 5:end))

qds = [2 2 2 3 3 3 4 4 4 6 6 6 8 8 8 10 10 10 12 12 12 14 14 14 16 16 16]

Y = tsne(fujimaki_19_mat')

gscatter(Y(:, 1), Y(:, 2), qds, [], [], 50)

% saving to ggplot

tsne_trans_t = table(Y(:, 1), Y(:, 2), qds', 'VariableNames', {'tsne_1', 'tsne_2', 'qds'})


%% generating data for metabolic atlas subsystem plots 

% lipid catabolism: representative subsystem - Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)
lipid_catabolism_rep_subsystem = "Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)"

idx_rep_subsystem = strcmp(lipid_catabolism_rep_subsystem, subsystems_m);

corr_r_lipid_catabolism_rep_subsystem = corr_r(idx_rep_subsystem)
rxn_name_rep_subsystem = model.rxns(idx_rep_subsystem)

lipid_catabolism_rep_subsystem_t = table(string(rxn_name_rep_subsystem), corr_r_lipid_catabolism_rep_subsystem)

lipid_catabolism_rep_subsystem_t.Properties.VariableNames = {'reaction', 'value'}


% lipid anabolism: representative subsystem - Fatty acid biosynthesis (odd-chain)
lipid_anabolism_rep_subsystem = "Fatty acid biosynthesis (odd-chain)"

idx_rep_subsystem_a = strcmp(lipid_anabolism_rep_subsystem, subsystems_m);

corr_r_lipid_anabolism_rep_subsystem = corr_r(idx_rep_subsystem_a)
rxn_name_rep_subsystem = model.rxns(idx_rep_subsystem_a)

lipid_anabolism_rep_subsystem_t = table(string(rxn_name_rep_subsystem), corr_r_lipid_anabolism_rep_subsystem)

lipid_anabolism_rep_subsystem_t.Properties.VariableNames = {'reaction', 'value'}


