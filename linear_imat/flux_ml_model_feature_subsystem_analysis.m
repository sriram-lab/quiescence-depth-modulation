%% Analyzes metabolic features of flux-based QDS predictor

%% Initializing Cobra + Model

clear;

% initializes cobra toolbox
initCobraToolbox;
changeCobraSolver('gurobi');

%% Analysis of linear regression model

% Loads model and fluxes
load QDS_model_3.mat

%% subsystem enrichment in flux-based model:

%% up with quiescence deepening:

% enrichment analysis:

% find subsystems in subset of reactions and converts to 1x1 string

% up in quiescence deepening
subsystems_c = (coef > 0)

subsystems_subset = string()

for i = 1:length(subsystems_c)
i
    subsystems_1 = (subsystems_c{i});
    subsystems_subset(i) = string(subsystems_1(1, 1));

end

% find total subsystems in all reactions and converts to 1x1 string
subsystems_c_m = model.subSystems
subsystems_m = string()

for i = 1:length(subsystems_c_m)
i
    subsystems_1_m = (subsystems_c_m{i});
    subsystems_m(i) = string(subsystems_1_m(1, 1));

end

subsystem_types = unique(subsystems_subset)
fold_enrichment = []
pval = []

for i = 1:length(subsystem_types)

    subsystem_types(i)

    num_sample = sum(strcmp(subsystems_subset, subsystem_types(i)))

    num_overall = sum(strcmp(subsystems_m, subsystem_types(i)))

    fold_enrichment(i) = (num_sample ./ length(subsystems_subset)) ./ (num_overall ./ length(model.subSystems))

    % p-value calculation
    pval(i) = 1 - hygecdf(num_sample, length(subsystems_m), num_overall, length(subsystems_subset))

end

% benjamini-hochberg FDR p-value correction
[h, p_crit, ~, qval] = fdr_bh(pval, 0.05)


pathway_analysis_table = table(subsystem_types', fold_enrichment', qval', 'VariableNames', {'Subsystem', 'Fold_enrichment', 'Qval'})

% sort by fold enrichment
pathway_analysis_table = sortrows(pathway_analysis_table, 2, 'descend')

% only keep q-value < 0.05
pathway_analysis_table = pathway_analysis_table(pathway_analysis_table.Qval < 0.05, :)

%% down with quiescence deepening:

% enrichment analysis:

% find subsystems in subset of reactions and converts to 1x1 string

% down in quiescence deepening
subsystems_c = (coef < 0)

subsystems_subset = string()

for i = 1:length(subsystems_c)

    subsystems_1 = (subsystems_c{i});
    subsystems_subset(i) = string(subsystems_1(1, 1));

end

% find total subsystems in all reactions and converts to 1x1 string
subsystems_c_m = model.subSystems
subsystems_m = string()

for i = 1:length(subsystems_c_m)

    subsystems_1_m = (subsystems_c_m{i});
    subsystems_m(i) = string(subsystems_1_m(1, 1));

end

subsystem_types = unique(subsystems_subset)
fold_enrichment = []
pval = []

for i = 1:length(subsystem_types)

    subsystem_types(i)

    num_sample = sum(strcmp(subsystems_subset, subsystem_types(i)))

    num_overall = sum(strcmp(subsystems_m, subsystem_types(i)))

    fold_enrichment(i) = (num_sample ./ length(subsystems_subset)) ./ (num_overall ./ length(model.subSystems))

    % p-value calculation
    pval(i) = 1 - hygecdf(num_sample, length(subsystems_m), num_overall, length(subsystems_subset))

end

% benjamini-hochberg FDR p-value correction
[h, p_crit, ~, qval] = fdr_bh(pval, 0.05)

pathway_analysis_table = table(subsystem_types', fold_enrichment', qval', 'VariableNames', {'Subsystem', 'Fold_enrichment', 'Qval'})

% sort by fold enrichment
pathway_analysis_table = sortrows(pathway_analysis_table, 2, 'descend')

% only keep q-value < 0.05
pathway_analysis_table = pathway_analysis_table(pathway_analysis_table.Qval < 0.05, :)




