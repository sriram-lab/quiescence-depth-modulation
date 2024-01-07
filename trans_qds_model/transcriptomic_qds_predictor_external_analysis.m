%% Applies multi-cell-type transcriptomics predictor to external datasets

%% Loads QDS model

clear
load multi_cell_type_transcriptomics_quiescence_depth_predictor.mat

%% johnson_18

johnson_18 = readtable("johnson_18_GSE117444_prolif_qui_count.csv")

% counts filter
counts_filter = sum(johnson_18{:, 2:end} >= 10, 2) >= 6 * 1;

johnson_18 = johnson_18(counts_filter, :);

% normalization: scaling each column by size factor
counts = []; normCounts = [];
counts = johnson_18{:, 2:end};

pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

normCounts = bsxfun(@rdivide,counts,sizeFactors);

johnson_18{:, 2:end} = normCounts;

johnson_18.gene = upper(johnson_18.gene);

johnson_18 = rmmissing(johnson_18)

% intersect
johnson_18_coef = innerjoin(coef_table, johnson_18);

coef_johnson = johnson_18_coef.coef;

johnson_18_mat = johnson_18_coef{:, 3:end};

% z-score
johnson_18_mat = normalize(johnson_18_mat, 2);

qds_pred = johnson_18_mat' * coef_johnson + coef0

[~, p] = ttest2(qds_pred(1:3), qds_pred(4:6))


%% clamer_18 - fpkm

clamer_18 = readtable("clamer_18_GSE112295_MCF7_fpkm_table.txt");

counts = []; normCounts = [];
counts = clamer_18{:, 4:end};

% normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

normCounts = bsxfun(@rdivide,counts,sizeFactors);

clamer_18{:, 4:end} = normCounts;

clamer_18 = rmmissing(clamer_18);

clamer_18.Properties.VariableNames{1} = 'gene'

% intersect
clamer_18_coef = innerjoin(coef_table, clamer_18);

coef_clamer = clamer_18_coef.coef;

clamer_18_mat = clamer_18_coef{:, 5:end};

% z-score
clamer_18_mat = normalize(clamer_18_mat, 2);
clamer_18_mat(isnan(clamer_18_mat)) = 0;

qds_pred = clamer_18_mat' * coef_clamer + coef0

% all: translatomics + mP control + transcriptomics
all_egf_idx = [9 10 15 16 3 4 21 22]
all_serum_starv_idx = [11 12 17 18 5 6 23 24]

mean(qds_pred(all_egf_idx))
mean(qds_pred(all_serum_starv_idx))

[~, p] = ttest2(qds_pred(all_egf_idx), qds_pred(all_serum_starv_idx))


%% min_18

min_18 = readtable("min_18_tpm_zscore.csv");
min_18 = rmmissing(min_18);

% intersect
min_18_coef = innerjoin(coef_table, min_18);

coef_min_18 = min_18_coef.coef;

min_18_mat = min_18_coef{:, 3:end};

qds_pred = min_18_mat' * coef_min_18 + coef0

ctrl_idx = [21 22 31 32]
p21_high_idx = [1 2 4 5 8 11 12 13 14 15]
p21_low_idx = [1 2 4 5 8 11 12 13 14 15]
serum_starv_idx = [23 24 33 34]
mek_inhib_idx = [25 26 35 36]
cdk46_inhib_idx = [27 28 37 38]
contact_inhib_idx = [29 30 39 40]

mean(qds_pred(ctrl_idx))

% serum starvation
[~, p] = ttest2(qds_pred(ctrl_idx), qds_pred(serum_starv_idx))

% spontaneous quiescence
[~, p] = ttest2(qds_pred(ctrl_idx), qds_pred(p21_high_idx))

% mek inhibition quiescence
[~, p] = ttest2(qds_pred(ctrl_idx), qds_pred(mek_inhib_idx))

% cdk46 inhibition quiescence
[~, p] = ttest2(qds_pred(ctrl_idx), qds_pred(cdk46_inhib_idx))

% contact inhibition quiescence
[~, p] = ttest2(qds_pred(ctrl_idx), qds_pred(contact_inhib_idx))


%% sharma_21

sharma_21 = readtable('sharma_21_prolif.csv');

sharma_21.Properties.VariableNames{1} = 'gene'

% counts filter: 
counts_filter = sum((sharma_21{:, 2:end} >= 10), 2) >= 16 * 1

sharma_21 = sharma_21(counts_filter, :)

counts = sharma_21{:, 2:end};

% normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

normCounts = bsxfun(@rdivide,counts,sizeFactors);

sharma_21{:, 2:end} = normCounts;

sharma_21 = rmmissing(sharma_21);

% intersect
sharma_21_coef = innerjoin(coef_table, sharma_21);

coef_sharma_21 = sharma_21_coef.coef;

sharma_21_mat = sharma_21_coef{:, 3:end};

% z-score
sharma_21_mat = normalize(sharma_21_mat, 2)

qds_pred = sharma_21_mat' * coef_sharma_21 + coef0

quiesc_idx = [8 10 12 14 1 3 6]
prolif_idx = [9 11 13 15 16 2 4 5 7]

mean(qds_pred(quiesc_idx))
mean(qds_pred(prolif_idx))

[~, p] = ttest2(qds_pred(quiesc_idx), qds_pred(prolif_idx))


%% kovatcheva_17 

kovatcheva_17 = readtable("kovatcheva_17_tpm.csv");

kovatcheva_17 = rmmissing(kovatcheva_17);

% intersect
kovatcheva_17_coef = innerjoin(coef_table, kovatcheva_17);

coef_kovatcheva_17 = kovatcheva_17_coef.coef;

kovatcheva_17_mat = kovatcheva_17_coef{:, 3:end};

kovatcheva_17_mat = normalize(kovatcheva_17_mat, 2)

qds_pred = kovatcheva_17_mat' * coef_kovatcheva_17 + coef0

[~, p] = ttest2(qds_pred(1:3), qds_pred(4:6))


%% kowalczyk_15

kowalczyk_15 = readtable("kowalczyk_15_bl6_log2_tpm_zscore.csv");

kowalczyk_15 = rmmissing(kowalczyk_15);

% intersect
kowalczyk_15_coef = innerjoin(coef_table, kowalczyk_15);

coef_kowalczyk_15 = kowalczyk_15_coef.coef;

kowalczyk_15_mat = kowalczyk_15_coef{:, 3:end};

qds_pred = kowalczyk_15_mat' * coef_kowalczyk_15 + coef0

pat = '_' + digitsPattern

grouping = kowalczyk_15.Properties.VariableNames(1:1434)
grouping = extractBefore(grouping, pat)

[~, p] = ttest2(qds_pred(strcmp('young_LT_HSC', grouping)), qds_pred(strcmp('young_ST_HSC', grouping))  )
[~, p] = ttest2(qds_pred(strcmp('young_LT_HSC', grouping)), qds_pred(strcmp('young_MPP', grouping))  )
[~, p] = ttest2(qds_pred(strcmp('young_ST_HSC', grouping)), qds_pred(strcmp('young_MPP', grouping))  )

[~, p] = ttest2(qds_pred(strcmp('old_LT_HSC', grouping)), qds_pred(strcmp('old_ST_HSC', grouping))  )
[~, p] = ttest2(qds_pred(strcmp('old_LT_HSC', grouping)), qds_pred(strcmp('old_MPP', grouping))  )

[~, p] = ttest2(qds_pred(strcmp('old_LT_HSC', grouping)), qds_pred(strcmp('young_LT_HSC', grouping))  )
[~, p] = ttest2(qds_pred(strcmp('old_ST_HSC', grouping)), qds_pred(strcmp('young_ST_HSC', grouping))  )
[~, p] = ttest2(qds_pred(strcmp('old_MPP', grouping)), qds_pred(strcmp('young_MPP', grouping))  )


%% GSE104651_apostolopoulou-svz-rna-seq.xlsx

apostolopoulou = readtable("GSE104651_apostolopoulou-svz-rna-seq.xlsx");

apostolopoulou = apostolopoulou(:, [2 4:end])

apostolopoulou.Properties.VariableNames{1} = 'gene'

apostolopoulou.gene = upper(apostolopoulou.gene)

apostolopoulou = rmmissing(apostolopoulou)

apostolopoulou_coef = innerjoin(coef_table, apostolopoulou);

coef_apostolopoulou = apostolopoulou_coef.coef;

apostolopoulou_mat = apostolopoulou_coef{:, 3:14};

apostolopoulou_mat = normalize(apostolopoulou_mat, 2)

apostolopoulou_mat(isnan(apostolopoulou_mat)) = 0

qds_pred = apostolopoulou_mat' * coef_apostolopoulou + coef0

age = [2 2 2 6 6 6 18 18 18 22 22 22]'

[r, p] = corr(age, qds_pred)

scatter(age, qds_pred)

% median
age_months_med = [2 6 18 22]

median_qds_pred = [median(qds_pred(1:3)), median(qds_pred(4:6)), median(qds_pred(7:9)), median(qds_pred(10:12))];

[r, p] = corr(age_months_med', median_qds_pred')

scatter(age_months_med, median_qds_pred)


%% maryanovich

maryanovich = readtable("GSE109546_HSC_Partek_Maria-aging_individual_counts_Maryanovich_18.xlsx")

% counts filter - all reads at least 10
counts_filter = sum(maryanovich{:, 2:end} >= 10, 2) >= 9 * 1;

maryanovich = maryanovich(counts_filter, :);

%% load maryanovich_specific_model
% model reconstructed with maryanovich-intersecting genes due to sparsity
% of model feature genes found in maryanovich

load maryanovich_model_v2_log_linear_optimized_elastic_net_model_3_datasets_normalized_qds.mat

counts = maryanovich{:, 2:end};

% MATLAB-recommended normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

normCounts = bsxfun(@rdivide,counts,sizeFactors);

maryanovich{:, 2:end} = normCounts;

maryanovich.Properties.VariableNames{1} = 'gene'

maryanovich.gene = upper(maryanovich.gene)

maryanovich = rmmissing(maryanovich);

% intersect
maryanovich_coef = innerjoin(coef_table, maryanovich);

coef_maryanovich = maryanovich_coef.coef;

maryanovich_mat = maryanovich_coef{:, 3:end};

maryanovich_mat = normalize(maryanovich_mat, 2);

maryanovich_mat(isnan(maryanovich_mat)) = 0;

qds_pred = maryanovich_mat' * coef_maryanovich + coef0

[~, p] = ttest2(qds_pred(1:3), qds_pred(4:6))


%% wang_21 embryogenesis

clear
load multi_cell_type_transcriptomics_quiescence_depth_predictor.mat

wang_21 = readtable("GSE136714_imputed_wang_21_embryogenesis_filtered_elastic_net.csv")
wang_21.gene = upper(wang_21.gene);

% omit missing genes
wang_21 = wang_21(~strcmp(wang_21.gene, ""), :)

counts = wang_21{:, 2:end};

% MATLAB-recommended normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

normCounts = bsxfun(@rdivide,counts,sizeFactors);
normCounts(1:10,:);

wang_21{:, 2:end} = normCounts;

% omit missing values
wang_21 = rmmissing(wang_21);

%% intersect
wang_21_coef = innerjoin(coef_table, wang_21);

coef_wang_21 = wang_21_coef.coef;

wang_21_mat = wang_21_coef{:, 3:end};

wang_21_mat = normalize(wang_21_mat, 2);

wang_21_mat(isnan(wang_21_mat)) = 0;

qds_pred = wang_21_mat' * coef_wang_21 + coef0

grouping = strings()
grouping(2:27) = "Zygote"
grouping(28:75) = "2 Cell"
grouping(76:169) = "4 Cell"
grouping(170:217) = "Late 4 Cell"
grouping(218:281) = "8 Cell"
grouping(282:329) = "16 Cell"
grouping(330:402) = "32 Cell"

grouping = grouping(2:end)

p = anova1(qds_pred, grouping)

[~, pval] = ttest2(qds_pred(strcmp(grouping, "Zygote")), qds_pred(strcmp(grouping, "2 Cell")))
[~, pval] = ttest2(qds_pred(strcmp(grouping, "2 Cell")), qds_pred(strcmp(grouping, "4 Cell")))
[~, pval] = ttest2(qds_pred(strcmp(grouping, "Late 4 Cell")), qds_pred(strcmp(grouping, "8 Cell")))
[~, pval] = ttest2(qds_pred(strcmp(grouping, "8 Cell")), qds_pred(strcmp(grouping, "16 Cell")))
[~, pval] = ttest2(qds_pred(strcmp(grouping, "16 Cell")), qds_pred(strcmp(grouping, "32 Cell")))


