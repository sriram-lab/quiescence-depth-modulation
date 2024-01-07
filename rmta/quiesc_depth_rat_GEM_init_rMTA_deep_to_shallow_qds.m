%% Predicting modulators of quiescence depth using INIT + rMTA

% Determines gene knockouts that move a cell from deep to shallow
% quiescence using the Robust Metabolic Transformation Algorithm

% NOTE: code based in part on COBRA rMTA tutorial:
% https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialRMTA.html

%% Initializing Cobra + Model

clear;

% initializes cobra toolbox
initCobraToolbox(false);
changeCobraSolver('gurobi');

%% Loads quiescence deepening transcriptomics
fujimaki_19 = readtable('fujimaki_19.txt');

%% loads Rat-GEM metabolic model
load Rat-GEM.mat;
model = ratGEM;
model.b = model.b(:, 1)
model = generateRules(model)
model.genes = upper(model.genes)

%% isolates genes in model

fujimaki_19.Gene = upper(fujimaki_19.Gene)
fujimaki_19 = fujimaki_19(ismember(fujimaki_19.Gene, model.genes), :)


%% normalizes data

% source state: deep quiescence (QDS of 12, 14) [22:27]
% target state: shallow quiescence (QDS of 2, 3) [4:9]

% isolate transcript levels
fujimaki_19_mat = table2array(fujimaki_19(:, 2:end));

% quantile normalize
norm_data = quantilenorm(fujimaki_19_mat);

% log2 normalize
norm_data = log2(norm_data + 1);

mean_all = mean(norm_data)
stdev_all = std(norm_data)

%% find high, medium, and low genes in source data

% source: deep quiescence (12-14 QDS)
source_mat = norm_data(:, 22:27)
    
alpha = 1.1; beta = 0.6;

%% run INIT to obtain metabolic model constrained by transcriptomics

deep_quiescence_model = run_init2(model, alpha, beta, source_mat, fujimaki_19.Gene, 34);

%% run CHRR to obtain flux samples

num_samples = 1E3

% runs CHRR
[v_ref, flux_all_array, rxn_inactive] = run_chrr2(model, deep_quiescence_model, num_samples);


%% calculate differentially expressed genes between source and target

% isolate transcript levels
fujimaki_19_mat = table2array(fujimaki_19(:, 2:end));

% quantile normalize
norm_data = quantilenorm(fujimaki_19_mat);

% re-generate source and target matrices w/o log2 normalization
source_mat = norm_data(:, 22:27)

target_mat = norm_data(:, 4:9)


gene_count_table = array2table([source_mat, target_mat])

gene_count_table.Properties.VariableNames = {'deep_quiescence_1', 'deep_quiescence_2', 'deep_quiescence_3', 'deep_quiescence_4', 'deep_quiescence_5', 'deep_quiescence_6',... 
                                              'shallow_quiescence_1', 'shallow_quiescence_2', 'shallow_quiescence_3', 'shallow_quiescence_4', 'shallow_quiescence_5', 'shallow_quiescence_6'}



gene_count_table.gene = string(fujimaki_19.Gene)

diff_table = rnaseqde(gene_count_table, ["shallow_quiescence_1", "shallow_quiescence_2", "shallow_quiescence_3", "shallow_quiescence_4", "shallow_quiescence_5", "shallow_quiescence_6"],...
                                         ["deep_quiescence_1", "deep_quiescence_2", "deep_quiescence_3", "deep_quiescence_4", "deep_quiescence_5", "deep_quiescence_6"],...
                     IDColumns= "gene");

figure
scatter(log2(mean([diff_table.Mean1, diff_table.Mean2], 2)), diff_table.Log2FoldChange, 3, diff_table.AdjustedPValue,'o')
colormap(flipud(cool(256)))
colorbar;
ylabel("log2(Fold change)")
xlabel("log2(Mean of normalized counts)")
title("Fold change by FDR")

% change column names so they work with diffexprs2rxnFBS
diff_table.logFC = diff_table.Log2FoldChange
diff_table.pval = diff_table.AdjustedPValue


%% find differential reactions

log_fc_thresh = 0.0; 
p_val_thresh = 5E-2; 

rxnFBS = diffexprs2rxnFBS(model, diff_table, v_ref, 'logFC', log_fc_thresh, 'pval', p_val_thresh);

% sets inactive reactions to zero
rxnFBS(rxn_inactive) = 0;

%% run rMTA

% Define alpha values to calculate rMTA
alpha_values = [0.66];

% Calculate epsilon, different for each reaction and with a minimum required change of 1e-3 (default)
epsilon = calculateEPSILON(flux_all_array, rxnFBS, 'minimum', 1E-3);

tic

[TSscore, deletedGenes, Vres] = rMTA(model, rxnFBS, v_ref, alpha_values, epsilon, ...
    'printLevel', 1, 'numWorkers', 36);

TIME.rMTA = toc  


%% table with results

rmta_results = table(deletedGenes, TSscore.bTS, TSscore.wTS, TSscore.mTS, TSscore.rTS)

% sort by rMTS
rmta_results = sortrows(rmta_results, 5, 'descend')

% omit missing values
rmta_results = rmmissing(rmta_results)

