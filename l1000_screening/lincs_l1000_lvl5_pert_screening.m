%% L1000 QDS Perturbation Screening

% Predicts effects of small molecule and genetic perturbations on
% quiescence depth using level 5 L1000 transcriptomics data

%% L1000 data loading

setup_env

gct_file_location = fullfile(cmapmpath, 'resources1', 'GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx');

lincs_pert_screen_level_5 = cmapm.Pipeline.parse_gctx(gct_file_location);


%% adds condition metadata

% conditions metadata for all treatment types
conditions_table = readtable("GSE70138_Broad_LINCS_sig_info.txt")

idx = [];

[~, idx] = ismember(lincs_pert_screen_level_5.cid(:, 1), string(conditions_table.sig_id))

lincs_pert_screen_level_5.cid(:, 2) = conditions_table.sig_id(idx);
lincs_pert_screen_level_5.cid(:, 3) = conditions_table.pert_id(idx);
lincs_pert_screen_level_5.cid(:, 4) = conditions_table.pert_iname(idx);
lincs_pert_screen_level_5.cid(:, 5) = conditions_table.pert_type(idx);
lincs_pert_screen_level_5.cid(:, 6) = conditions_table.cell_id(idx);
lincs_pert_screen_level_5.cid(:, 7) = conditions_table.pert_idose(idx);
lincs_pert_screen_level_5.cid(:, 8) = conditions_table.pert_itime(idx);
lincs_pert_screen_level_5.cid(:, 9) = conditions_table.distil_id(idx);

% confirm conditions are equivalent
isequal(lincs_pert_screen_level_5.cid(:, 1), lincs_pert_screen_level_5.cid(:, 2))

%% adds gene metadata

% gene metadata 
gene_table = readtable("GSE92742_Broad_LINCS_gene_info.txt")

idx = [];

[~, idx] = ismember(lincs_pert_screen_level_5.rid, cellstr(string(gene_table.pr_gene_id)))

lincs_pert_screen_level_5.rid(:, 2) = cellstr(string(gene_table.pr_gene_id(idx)));
lincs_pert_screen_level_5.rid(:, 3) = gene_table.pr_gene_symbol(idx);
lincs_pert_screen_level_5.rid(:, 4) = cellstr(string(gene_table.pr_is_lm(idx)));
lincs_pert_screen_level_5.rid(:, 5) = cellstr(string(gene_table.pr_is_bing(idx)));

% confirm conditions are equivalent
isequal(lincs_pert_screen_level_5.rid(:, 1), lincs_pert_screen_level_5.rid(:, 2))

%% saves to .mat file

% save("lincs_pert_screen_level_5.mat", "lincs_pert_screen_level_5", "-v7.3")

%% gene knockout screening

% load l1000 data
load lincs_pert_screen_level_5.mat

%% load elastic net model

load l1000_landmark_model_v2_log_linear_optimized_elastic_net_model_3_datasets_normalized_qds.mat

transcriptomics.Properties.VariableNames(1) = {'genes'}
    
transcriptomics.genes = string(transcriptomics.genes)

clear genes


%% rows: isolates genes that are in elastic net model training set

genes = strings();

lincs_gene_screen_mat = lincs_pert_screen_level_5.mat(ismember(lincs_pert_screen_level_5.rid(:, 3), transcriptomics.genes), :);

genes = lincs_pert_screen_level_5.rid(ismember(lincs_pert_screen_level_5.rid(:, 3), transcriptomics.genes), 3);


%% isolate landmark genes (directly measured)

lincs_gene_screen_mat = lincs_pert_screen_level_5.mat(ismember(lincs_pert_screen_level_5.rid(:, 3), transcriptomics.genes) & strcmp(lincs_pert_screen_level_5.rid(:, 4), "1"), :);

genes = lincs_pert_screen_level_5.rid(ismember(lincs_pert_screen_level_5.rid(:, 3), transcriptomics.genes) & strcmp(lincs_pert_screen_level_5.rid(:, 4), "1"), 3);


%% columns: isolates CRISPR perturbagens 

% transpose so it aligns with columns (unnecessary if cid already aligned
% with columns)
lincs_pert_screen_level_5.cid = lincs_pert_screen_level_5.cid'

intervention_type = strings(); cell_line = strings();

intervention_type = string(lincs_pert_screen_level_5.cid(5, :));

cell_line = string(lincs_pert_screen_level_5.cid(6, :));


%% types of cells available:
unique(cell_line(strcmp(intervention_type, "trt_xpr")))


%% experimental conditions

lincs_gene_screen_mat_exp = lincs_gene_screen_mat(:, strcmp(intervention_type, "trt_xpr"));

conditions_exp = lincs_pert_screen_level_5.cid(:, strcmp(intervention_type, "trt_xpr"));


%% control conditions (3 types: untreated, vector, and vehicle)

% vector: best control for genetic perturbations
lincs_gene_screen_mat_ctl_vector = lincs_gene_screen_mat(:, strcmp(intervention_type, "ctl_vector"));

conditions_ctl_vector = lincs_pert_screen_level_5.cid(:, strcmp(intervention_type, "ctl_vector"));


%% gene KO analysis

gene_ko_unique = unique(conditions_exp(4, :));
gene_ko_name = strings()
pval = [];


%% predict QDS for control: vector

% calculates qds for control conditions

% creates table of gene name + transcriptomics for control conditions
ctrl_table = array2table(lincs_gene_screen_mat_ctl_vector);
ctrl_table.genes = genes;

% regression model table
regression_model_table = table(transcriptomics{:, 1}, coef);

regression_model_table.Properties.VariableNames = {'genes', 'coef'};

% merges gene_table with regression model table (genes column + weights)
merged_table_ctrl = innerjoin(ctrl_table, regression_model_table, 'Keys', {'genes'});

% creates new coefficient vector genes
coef_i_ctrl = merged_table_ctrl.coef;

transcriptomics_mat_ctrl = merged_table_ctrl{:, 1:width(lincs_gene_screen_mat_ctl_vector)};

width(lincs_gene_screen_mat_ctl_vector)

% qds prediction array
qds_pred_ctrl = [];

qds_pred_ctrl = transcriptomics_mat_ctrl' * coef_i_ctrl + coef0;

mean(qds_pred_ctrl)


%% predicts qds for each gene ko

for i = 1:length(gene_ko_unique)

  
    i 
    sig_ko = [];
    idx_ko = [];

    idx_ko = strcmp(gene_ko_unique{i}, conditions_exp(4, :));

    % variable to verify conditions are correct
    ko_ver = conditions_exp(:, idx_ko);

    % isolates transcriptomics for ko
    sig_ko = lincs_gene_screen_mat_exp(:, idx_ko);

    % creates table of gene name + transcriptomics for ko
    ko_table = array2table(sig_ko);
    ko_table.genes = genes;

    % regression model table
    regression_model_table = table(transcriptomics{:, 1}, coef);

    regression_model_table.Properties.VariableNames = {'genes', 'coef'};
    
    % merges ko table with regression model table (genes column + weights)
    merged_table = innerjoin(ko_table, regression_model_table, 'Keys', {'genes'});

    % creates new coefficient vector 
    coef_i = merged_table.coef;

    transcriptomics_mat = merged_table{:, 1:width(sig_ko)};

    % qds prediction array
    qds_pred = [];

    qds_pred = transcriptomics_mat' * coef_i + coef0;

    gene_ko_name(i) = gene_ko_unique{i};
    qds_pred_mean(i) = mean(qds_pred);
    qds_pred_over_ctrl_mean(i) = mean(qds_pred) / mean(qds_pred_ctrl);

    % two-sample t-test: ko vs control
    [h, pval(i)] = ttest2(qds_pred, qds_pred_ctrl);

end

%% aggregates results

% benjamini-hochberg p-value correction
[h, p_crit, ~, p_adj] = fdr_bh(pval)

gene_screen_lincs_2_table = table(gene_ko_unique', qds_pred_mean', qds_pred_over_ctrl_mean', pval', p_adj', 'VariableNames', {'gene_ko_name', 'qds_pred_mean', 'qds_pred_over_ctrl_mean', 'p_value', 'q_value'})

gene_screen_lincs_2_table = sortrows(gene_screen_lincs_2_table, 3)


%% Small molecule screening

clear

% load l1000 data
load lincs_pert_screen_level_5.mat

%% loads elastic net model

% model built using landmark genes
load l1000_landmark_model_v2_log_linear_optimized_elastic_net_model_3_datasets_normalized_qds.mat

transcriptomics.Properties.VariableNames(1) = {'genes'}
    
transcriptomics.genes = string(transcriptomics.genes)

clear genes


%% rows: isolates genes that are in elastic net model

genes = strings();

lincs_drug_screen_mat = lincs_pert_screen_level_5.mat(ismember(lincs_pert_screen_level_5.rid(:, 3), transcriptomics.genes), :);

genes = lincs_pert_screen_level_5.rid(ismember(lincs_pert_screen_level_5.rid(:, 3), transcriptomics.genes), 3);

%%
dose = []; time = []; cell_line = strings(); pert_type = string();

lincs_pert_screen_level_5.cid = lincs_pert_screen_level_5.cid'

dose = (lincs_pert_screen_level_5.cid(7, :));
time = (lincs_pert_screen_level_5.cid(8, :));
cell_line = string(lincs_pert_screen_level_5.cid(6, :));
pert_type = string(lincs_pert_screen_level_5.cid(5, :));


%% isolate experimental conditions
lincs_2_sig_experiment = lincs_drug_screen_mat(:, strcmp(pert_type, "trt_cp") & strcmp(dose, "10.0 um") & strcmp(time, "24 h"));

conditions_experiment = lincs_pert_screen_level_5.cid(:, strcmp(pert_type, "trt_cp") & strcmp(dose, "10.0 um") & strcmp(time, "24 h"));

%% isolate control conditions

% vehicle control is most appropriate since treatment is compound (not
% genetic perturbation)
lincs_2_sig_control = lincs_drug_screen_mat(:, strcmp(pert_type, "ctl_vehicle") & strcmp(time, "24 h"));

conditions_control = lincs_pert_screen_level_5.cid(:, strcmp(pert_type, "ctl_vehicle") & strcmp(time, "24 h"));


%% drug analysis
% for each drug, calculate qds for treatment and control
% then compare qds with paired t-test and find ratio of q_t / c_ctr

drugs_unique = unique(conditions_experiment(4, :));

drug_name = strings()
qds_pred_mean = zeros([length(drugs_unique) 1]);
qds_pred_over_ctrl_mean = zeros([length(drugs_unique) 1]);
pval = zeros([length(drugs_unique) 1]);


%% calculates qds for control conditions

% creates table of gene name + transcriptomics for control conditions
% (DMSO, same 24 hour time point)
ctrl_table = array2table(lincs_2_sig_control);
ctrl_table.genes = genes;

% regression model table
regression_model_table = table(transcriptomics{:, 1}, coef);

regression_model_table.Properties.VariableNames = {'genes', 'coef'};

% merges drug_table with regression model table (genes column + weights)
merged_table_ctrl = innerjoin(ctrl_table, regression_model_table, 'Keys', {'genes'});

% creates new coefficient vector with repeated coefficients for
% repeated genes
coef_i_ctrl = merged_table_ctrl.coef;

transcriptomics_mat_ctrl = merged_table_ctrl{:, 1:width(lincs_2_sig_control)};

width(lincs_2_sig_control)

% qds prediction array
qds_pred_ctrl = [];

qds_pred_ctrl = transcriptomics_mat_ctrl' * coef_i_ctrl + coef0;

mean(qds_pred_ctrl)


%% predicts qds for each small molecule

for i = 1:length(drugs_unique)

    sig_drug = [];
    idx_drug = [];

    % drug
    idx_drug = strcmp(drugs_unique{i}, conditions_experiment(4, :));

    % vector to verify same drug is selected
    drug_ver = conditions_experiment(7, idx_drug);

    % isolates transcriptomics for drug
    sig_drug = lincs_2_sig_experiment(:, idx_drug);

    % creates table of gene name + transcriptomics for drug
    drug_table = array2table(sig_drug);
    drug_table.genes = genes;

    % regression model table
    regression_model_table = table(transcriptomics{:, 1}, coef);

    regression_model_table.Properties.VariableNames = {'genes', 'coef'};
    
    % merges drug_table with regression model table (genes column + weights)
    merged_table = innerjoin(drug_table, regression_model_table, 'Keys', {'genes'});

    % creates new coefficient vector with repeated coefficients for
    % repeated genes
    coef_i = merged_table.coef;

    transcriptomics_mat = merged_table{:, 1:width(sig_drug)};

    % qds prediction array
    qds_pred = [];

    qds_pred = transcriptomics_mat' * coef_i + coef0;
   
    drug_name(i) = drugs_unique{i};
    qds_pred_mean(i) = mean(qds_pred);
    qds_pred_over_ctrl_mean(i) = mean(qds_pred) / mean(qds_pred_ctrl);

    % two-sample t-test: drug vs control
    [h, pval(i)] = ttest2(qds_pred, qds_pred_ctrl);
    
end


%% aggregates results

% benjamini-hochberg p-value correction
[h, p_crit, ~, qval] = fdr_bh(pval)

drug_screen_lincs_2_table = table(drug_name', qds_pred_mean, qds_pred_over_ctrl_mean, pval, qval, 'VariableNames', {'drug_name', 'qds_pred_mean', 'qds_pred_over_ctrl_mean', 'p_value', 'q_value'})

drug_screen_lincs_2_table = sortrows(drug_screen_lincs_2_table, 3, 'descend')


