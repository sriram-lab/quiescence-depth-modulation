%% External validation of flux-based QDS predictor 

% Predicts QDS for several external datasets using flux-based QDS predictor

%% Initializes Cobra 

clear;

% initializes cobra toolbox
initCobraToolbox;
changeCobraSolver('gurobi');

% loads Rat-GEM, optimized regression model, optimized params
load QDS_model_3.mat


%% testing model on fluxes generated through INIT + CHRR

load init_chrr_alpha_1.1_beta_0.6_all_qds.MAT

flux_init_chrr = cell2mat(flux_avg_array);

%% results

qds_pred_init_chrr = flux_init_chrr' * coef + coef0

qds_avg = [2 3 4 6 8 10 12 14 16]

[r, p] = corr(qds_avg', qds_pred_init_chrr)

scatter(qds_avg', qds_pred_init_chrr)

init_chrr_val_table = table(qds_avg', qds_pred_init_chrr, 'VariableNames', {'qds_actual', 'qds_pred'})


%% External validation:

%% sharma_21: proliferation vs quiesence 

sharma_21 = readtable('sharma_21_prolif.csv');

% counts filter: all reads ≥ 10 count
counts_filter = sum((sharma_21{:, 2:end} >= 10), 2) >= 16 * 1

sharma_21 = sharma_21(counts_filter, :)

counts = sharma_21{:, 2:end};

% MATLAB-recommended normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

normCounts = bsxfun(@rdivide,counts,sizeFactors);

sharma_21{:, 2:end} = normCounts;

transcriptomics = sharma_21;

% Isolating Conditions

thresh_up = 1;

thresh_down = 1 / thresh_up;

EGF1_up = intersect(transcriptomics.Var1(transcriptomics.c2 ./ transcriptomics.c1 > thresh_up), model.genes);
EGF1_down = intersect(transcriptomics.Var1(transcriptomics.c2 ./ transcriptomics.c1 < thresh_down), model.genes);

EGF2_up = intersect(transcriptomics.Var1(transcriptomics.c6 ./ transcriptomics.c5 > thresh_up), model.genes);
EGF2_down = intersect(transcriptomics.Var1(transcriptomics.c6 ./ transcriptomics.c5 < thresh_down), model.genes);

EGF3_up = intersect(transcriptomics.Var1(transcriptomics.c10 ./ transcriptomics.c9 > thresh_up), model.genes);
EGF3_down = intersect(transcriptomics.Var1(transcriptomics.c10 ./ transcriptomics.c9 < thresh_down), model.genes);

TPA1_up = intersect(transcriptomics.Var1(transcriptomics.c13 ./ transcriptomics.c5 > thresh_up), model.genes);
TPA1_down = intersect(transcriptomics.Var1(transcriptomics.c13 ./ transcriptomics.c5 < thresh_down), model.genes);

TPA2_up = intersect(transcriptomics.Var1(transcriptomics.c16 ./ transcriptomics.c15 > thresh_up), model.genes);
TPA2_down = intersect(transcriptomics.Var1(transcriptomics.c16 ./ transcriptomics.c15 < thresh_down), model.genes);

H89_EGF1_up = intersect(transcriptomics.Var1(transcriptomics.c8 ./ transcriptomics.c7 > thresh_up), model.genes)
H89_EGF1_down = intersect(transcriptomics.Var1(transcriptomics.c8 ./ transcriptomics.c7 < thresh_down), model.genes)

H89_EGF2_up = intersect(transcriptomics.Var1(transcriptomics.c12 ./ transcriptomics.c11 > thresh_up), model.genes)
H89_EGF2_down = intersect(transcriptomics.Var1(transcriptomics.c12 ./ transcriptomics.c11 < thresh_down), model.genes)

H89_TPA1_up = intersect(transcriptomics.Var1(transcriptomics.c14 ./ transcriptomics.c7 > thresh_up), model.genes)
H89_TPA1_down = intersect(transcriptomics.Var1(transcriptomics.c14 ./ transcriptomics.c7 < thresh_down), model.genes)

H89_TPA2_up = intersect(transcriptomics.Var1(transcriptomics.c18 ./ transcriptomics.c17 > thresh_up), model.genes)
H89_TPA2_down = intersect(transcriptomics.Var1(transcriptomics.c18 ./ transcriptomics.c17 < thresh_down), model.genes)


up = {}; down = {};

up = {EGF1_up, EGF2_up, EGF3_up, TPA1_up, TPA2_up, H89_EGF1_up, H89_EGF2_up, H89_TPA1_up, H89_TPA2_up};
down = {EGF1_down, EGF2_down, EGF3_down, TPA1_down, TPA2_down, H89_EGF1_down, H89_EGF2_down, H89_TPA1_down, H89_TPA2_down};

num_rxns = length(model.rxnNames);

flux_up_sharma_21 = [];
flux_dw_sharma_21 = [];

parfor i = 1:length(up)
  
    flux_up_sharma_21(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    
    flux_dw_sharma_21(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);
 
end


%% results

qds_up = []; qds_dw = [];

% elastic net
qds_up = flux_up_sharma_21 * coef + coef0
qds_dw = flux_dw_sharma_21 * coef + coef0
 
[h, pval] = ttest(qds_up, qds_dw)

mean(qds_up)
mean(qds_dw)

[h, pval] = ttest(qds_up(1:3), qds_dw(1:3))
[h, pval] = ttest(qds_up(1:5), qds_dw(1:5))
[h, pval] = ttest(qds_up(1:9), qds_dw(1:9))

qds_all = [qds_up; qds_dw]
grouping = strings()
grouping(1:9) = "Proliferation"
grouping(10:18) = "Quiescence"

sharma_21_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% johnson_18: proliferation vs quiesence (contact-inhibition-induced)

johnson_18 = readtable('johnson_18_prolif_qui_count.csv');

% counts filter - all reads at least 10, 100% of the time
counts_filter = sum(johnson_18{:, 2:end} >= 10, 2) >= 6 * 1

johnson_18 = johnson_18(counts_filter, :)
counts = []; normCounts = [];
counts = johnson_18{:, 2:end};

% normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

normCounts = bsxfun(@rdivide,counts,sizeFactors);

johnson_18{:, 2:end} = normCounts;

% isolates genes in model
johnson_18 = johnson_18(ismember(johnson_18.gene, model.genes), :);


% determines up / down genes
thresh_up = 1;
thresh_down = 1 / thresh_up;

up = {}; down = {};

up = {johnson_18.gene(johnson_18.P_rep1 ./ johnson_18.x7dCI_rep1 > thresh_up),...
      johnson_18.gene(johnson_18.P_rep2 ./ johnson_18.x7dCI_rep2 > thresh_up),...
      johnson_18.gene(johnson_18.P_rep3 ./ johnson_18.x7dCI_rep3 > thresh_up)
     };

down = {johnson_18.gene(johnson_18.P_rep1 ./ johnson_18.x7dCI_rep1 < thresh_down),...
      johnson_18.gene(johnson_18.P_rep2 ./ johnson_18.x7dCI_rep2 < thresh_down),...
      johnson_18.gene(johnson_18.P_rep3 ./ johnson_18.x7dCI_rep3 < thresh_down)
     };

flux_up_johnson_18 = [];
flux_dw_johnson_18 = [];

parfor i = 1:length(up)
  
    flux_up_johnson_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    
    flux_dw_johnson_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

end


%% results

qds_up = []; qds_dw = [];

% elastic net
qds_up = flux_up_johnson_18 * coef + coef0
qds_dw = flux_dw_johnson_18 * coef + coef0


[h, pval] = ttest(qds_up, qds_dw)

qds_all = [qds_up; qds_dw]
grouping = strings()
grouping(1:3) = "Proliferation"
grouping(4:6) = "Quiescence (Contact-Inhibition)"

johnson_18_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% clamer 18: proliferation vs quiesence vs control

clamer_18 = readtable('clamer_18_GSE112295_MCF7_fpkm_table.txt');

counts = []; normCounts = [];
counts = clamer_18{:, 4:end};

% MATLAB-recommended normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

normCounts = bsxfun(@rdivide,counts,sizeFactors);

clamer_18{:, 4:end} = normCounts;

% isolates genes in model
clamer_18 = clamer_18(ismember(clamer_18.gene_short_name, model.genes), :);


% determines up / down genes
thresh_up = 1;
thresh_down = 1 / thresh_up;


up = {}; down = {};

up = {clamer_18.gene_short_name(clamer_18.total_EGF_1 ./ clamer_18.total_STARV_1 > thresh_up),...
      clamer_18.gene_short_name(clamer_18.total_EGF_2 ./ clamer_18.total_STARV_2 > thresh_up),...
      clamer_18.gene_short_name(clamer_18.ribolace_EGF_1 ./ clamer_18.ribolace_STARV_1 > thresh_up),...
      clamer_18.gene_short_name(clamer_18.ribolace_EGF_2 ./ clamer_18.ribolace_STARV_2 > thresh_up),...
      clamer_18.gene_short_name(clamer_18.polysomal_EGF_1 ./ clamer_18.polysomal_STARV_1 > thresh_up),...
      clamer_18.gene_short_name(clamer_18.polysomal_EGF_2 ./ clamer_18.polysomal_STARV_2 > thresh_up),...
      clamer_18.gene_short_name(clamer_18.mP_EGF_1 ./ clamer_18.mP_STARV_1 > thresh_up),...
      clamer_18.gene_short_name(clamer_18.mP_EGF_2 ./ clamer_18.mP_STARV_2 > thresh_up)
     };

down = {clamer_18.gene_short_name(clamer_18.total_EGF_1 ./ clamer_18.total_STARV_1 < thresh_down),...
        clamer_18.gene_short_name(clamer_18.total_EGF_2 ./ clamer_18.total_STARV_2 < thresh_down),...
        clamer_18.gene_short_name(clamer_18.ribolace_EGF_1 ./ clamer_18.ribolace_STARV_1 < thresh_down),...
        clamer_18.gene_short_name(clamer_18.ribolace_EGF_2 ./ clamer_18.ribolace_STARV_2 < thresh_down),...
        clamer_18.gene_short_name(clamer_18.polysomal_EGF_1 ./ clamer_18.polysomal_STARV_1 < thresh_down),...
        clamer_18.gene_short_name(clamer_18.polysomal_EGF_2 ./ clamer_18.polysomal_STARV_2 < thresh_down),...
        clamer_18.gene_short_name(clamer_18.mP_EGF_1 ./ clamer_18.mP_STARV_1 < thresh_down),...
        clamer_18.gene_short_name(clamer_18.mP_EGF_2 ./ clamer_18.mP_STARV_2 < thresh_down)
     };


flux_up_clamer_18 = [];
flux_dw_clamer_18 = [];

parfor i = 1:length(up)
  

    i
    flux_up_clamer_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    
    flux_dw_clamer_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

    
end


%% results
qds_up = []; qds_dw = [];

% elastic net
qds_up = flux_up_clamer_18 * coef + coef0
qds_dw = flux_dw_clamer_18 * coef + coef0

[h, pval] = ttest(qds_up, qds_dw)

% exclude transcriptomics and mP control (just translatomics)
qds_up = flux_up_clamer_18(3:6, :) * coef + coef0
qds_dw = flux_dw_clamer_18(3:6, :) * coef + coef0

[h, pval] = ttest(qds_up, qds_dw)

qds_all = [qds_up; qds_dw]
grouping = strings()
grouping(1:4) = "Proliferation"
grouping(5:8) = "Quiescence"

clamer_18_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% min 18: proliferation vs quiesence (5 types) vs control

min_18 = readtable('min_18_GSE122927_ReadCount.csv');

% count filter
counts_filter = sum(min_18{:, 2:end} >= 10, 2) >= 40 * 1

min_18 = min_18(counts_filter, :)

counts = []; normCounts = [];
counts = min_18{:, 2:end};

% normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

normCounts = bsxfun(@rdivide,counts,sizeFactors);

min_18{:, 2:end} = normCounts;

% isolates genes in model
min_18 = min_18(ismember(min_18.gene, model.genes), :);

thresh_up = 1
thresh_down = 1 / thresh_up;


%% serum starvation

up = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.SerumStarvation_2E2_rep1 > thresh_up),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.SerumStarvation_2E2_rep2 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.SerumStarvation_3B6_rep1 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.SerumStarvation_3B6_rep2 > thresh_up),...
     };

down = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.SerumStarvation_2E2_rep1 < thresh_down),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.SerumStarvation_2E2_rep2 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.SerumStarvation_3B6_rep1 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.SerumStarvation_3B6_rep2 < thresh_down),...
     };


flux_up_min_18 = [];
flux_dw_min_18 = [];

parfor i = 1:length(up)
  

    i
    flux_up_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    
    flux_dw_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

    
end


qds_up = []; qds_dw = [];

% elastic net
qds_up = flux_up_min_18 * coef + coef0
qds_dw = flux_dw_min_18 * coef + coef0

[h, pval] = ttest(qds_up, qds_dw)

qds_all = [qds_up; qds_dw]
grouping = strings()
grouping(1:4) = "Proliferation"
grouping(5:8) = "Quiescence"

min_18_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% contact inhibition

thresh_up = 1
thresh_down = 1 / thresh_up;

up = {}; down = {};

up = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.ContactIn_2E2_rep1 > thresh_up),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.ContactIn_2E2_rep2 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.ContactIn_3B6_rep1 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.ContactIn_3B6_rep2 > thresh_up),...
     };

down = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.ContactIn_2E2_rep1 < thresh_down),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.ContactIn_2E2_rep2 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.ContactIn_3B6_rep1 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.ContactIn_3B6_rep2 < thresh_down),...
     };


flux_up_min_18 = [];
flux_dw_min_18 = [];

parfor i = 1:length(up)
  
    flux_up_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    
    flux_dw_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

    
end

qds_up = []; qds_dw = [];

% elastic net
qds_up = flux_up_min_18 * coef + coef0
qds_dw = flux_dw_min_18 * coef + coef0

[h, pval] = ttest(qds_up, qds_dw)


qds_all = [qds_up; qds_dw]
grouping = strings()
grouping(1:4) = "Proliferation"
grouping(5:8) = "Quiescence"

min_18_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% mek inhibition

thresh_up = 1
thresh_down = 1 / thresh_up;

up = {}; down = {};

up = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.Meki_2E2_rep1 > thresh_up),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.Meki_2E2_rep2 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.Meki_3B6_rep1 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.Meki_3B6_rep2 > thresh_up),...
     };

down = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.Meki_2E2_rep1 < thresh_down),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.Meki_2E2_rep2 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.Meki_3B6_rep1 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.Meki_3B6_rep2 < thresh_down),...
     };


flux_up_min_18 = [];
flux_dw_min_18 = [];

parfor i = 1:length(up)
  
    flux_up_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    
    flux_dw_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

    
end


qds_up = []; qds_dw = [];

% elastic net
qds_up = flux_up_min_18 * coef + coef0
qds_dw = flux_dw_min_18 * coef + coef0

[h, pval] = ttest(qds_up, qds_dw)


qds_all = [qds_up; qds_dw]
grouping = strings()
grouping(1:4) = "Proliferation"
grouping(5:8) = "Quiescence"

min_18_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% cdk46 inhibition

thresh_up = 1
thresh_down = 1 / thresh_up

up = {}; down = {};

up = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.CDK46i_2E2_rep1 > thresh_up),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.CDK46i_2E2_rep2 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.CDK46i_3B6_rep1 > thresh_up),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.CDK46i_3B6_rep2 > thresh_up),...
     };

down = {min_18.gene(min_18.Control_2E2_rep1 ./ min_18.CDK46i_2E2_rep1 < thresh_down),...
      min_18.gene(min_18.Control_2E2_rep2 ./ min_18.CDK46i_2E2_rep2 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep1 ./ min_18.CDK46i_3B6_rep1 < thresh_down),...
      min_18.gene(min_18.Control_3B6_rep2 ./ min_18.CDK46i_3B6_rep2 < thresh_down),...
     };

flux_up_min_18 = [];
flux_dw_min_18 = [];

parfor i = 1:length(up)
  

    i
    flux_up_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    
    flux_dw_min_18(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

    
end


qds_up = []; qds_dw = [];

% elastic net
qds_up = flux_up_min_18 * coef + coef0
qds_dw = flux_dw_min_18 * coef + coef0

[h, pval] = ttest(qds_up, qds_dw)


qds_all = [qds_up; qds_dw]
grouping = strings()
grouping(1:4) = "Proliferation"
grouping(5:8) = "Quiescence"

min_18_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% lee_17 - mmc2 multiomics analysis

lee_17_proteomics = readtable('mmc2.xlsx', 'Sheet', 'Proteins')
lee_17_proteomics.Gene = upper(lee_17_proteomics.Gene)

counts = []; normCounts = [];
counts = lee_17_proteomics{:, 4:9};


% normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

normCounts = bsxfun(@rdivide,counts,sizeFactors);

lee_17_proteomics{:, 4:9} = normCounts;

% isolates genes in model
lee_17_proteomics = lee_17_proteomics(ismember(lee_17_proteomics.Gene, model.genes), :);


% determines up / down genes
thresh_up = 1;
thresh_dw = 1 / thresh_up;


up = {}; down = {};

% comparison analysis: proliferative (8, 12, 16, and 20H) vs quiescent (0H)
up = {lee_17_proteomics.Gene(lee_17_proteomics.T_8H ./ lee_17_proteomics.T_0H > thresh_up), ...
      lee_17_proteomics.Gene(lee_17_proteomics.T_12H ./ lee_17_proteomics.T_0H > thresh_up), ...
      lee_17_proteomics.Gene(lee_17_proteomics.T_16H ./ lee_17_proteomics.T_0H > thresh_up), ...
      lee_17_proteomics.Gene(lee_17_proteomics.T_20H ./ lee_17_proteomics.T_0H > thresh_up)}

dw = {lee_17_proteomics.Gene(lee_17_proteomics.T_8H ./ lee_17_proteomics.T_0H < thresh_dw), ...
      lee_17_proteomics.Gene(lee_17_proteomics.T_12H ./ lee_17_proteomics.T_0H < thresh_dw), ...
      lee_17_proteomics.Gene(lee_17_proteomics.T_16H ./ lee_17_proteomics.T_0H < thresh_dw), ...
      lee_17_proteomics.Gene(lee_17_proteomics.T_20H ./ lee_17_proteomics.T_0H < thresh_dw)}

flux_lee_17_proteomics_up = [];
flux_lee_17_proteomics_dw = [];

parfor i = 1:length(up)
  
    i
    flux_lee_17_proteomics_up(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, dw{i}, kappa, rho, epsilon, 0, [], 1);
    
    flux_lee_17_proteomics_dw(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  dw{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

end

%% results

qds_lee_17_proteomics_up = flux_lee_17_proteomics_up * coef + coef0
qds_lee_17_proteomics_dw = flux_lee_17_proteomics_dw * coef + coef0

[h, p] = ttest(qds_lee_17_proteomics_up, qds_lee_17_proteomics_dw)


qds_all = [qds_lee_17_proteomics_up; qds_lee_17_proteomics_dw]
grouping = strings()
grouping(1:4) = "Proliferation"
grouping(5:8) = "Quiescence"

% saving results to ggplot
qds_lee_17_proteomics_t = table(qds_all, grouping', 'VariableNames', {'QDS', 'grouping'})


%% kovatcheva_17

kovatcheva_17 = readtable("kovatcheva_17_counts_with_gene_symbol.csv");

kovatcheva_17 = rmmissing(kovatcheva_17);

kovatcheva_17 = kovatcheva_17(ismember(kovatcheva_17.gene, model.genes), :)

kovatcheva_17.Properties.VariableNames(1:6) = {'serum_starv_1', 'serum_starv_2', 'serum_starv_3', 'prolif_1', 'prolif_2', 'prolif_3'}


counts = []; normCounts = [];
counts = kovatcheva_17{:, 1:6};


% normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

normCounts = bsxfun(@rdivide,counts,sizeFactors);

kovatcheva_17{:, 1:6} = normCounts;

%%
thresh_up = 1;

thresh_down = 1 / thresh_up;

up = {}; down = {};

up = {kovatcheva_17.gene(kovatcheva_17.prolif_1 ./ kovatcheva_17.serum_starv_1 > thresh_up),...
      kovatcheva_17.gene(kovatcheva_17.prolif_2 ./ kovatcheva_17.serum_starv_2 > thresh_up),...
      kovatcheva_17.gene(kovatcheva_17.prolif_3 ./ kovatcheva_17.serum_starv_3 > thresh_up)
     };

down = {kovatcheva_17.gene(kovatcheva_17.prolif_1 ./ kovatcheva_17.serum_starv_1 < thresh_down),...
      kovatcheva_17.gene(kovatcheva_17.prolif_2 ./ kovatcheva_17.serum_starv_2 < thresh_down),...
      kovatcheva_17.gene(kovatcheva_17.prolif_3 ./ kovatcheva_17.serum_starv_3 < thresh_down)
     };


flux_kovatcheva_17 = [];

parfor i = 1:length(up)
  

    flux_up_kovatcheva_17(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up{i}, down{i}, kappa, rho, epsilon, 0, [], 1);
    
    flux_dw_kovatcheva_17(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  down{i}, up{i}, kappa, rho, epsilon, 0, [], 1);

    
end

%% results

qds_pred_up = []; qds_pred_dw = [];

qds_pred_up = flux_up_kovatcheva_17 * coef + coef0

qds_pred_dw = flux_dw_kovatcheva_17 * coef + coef0

[h, p] = ttest(qds_pred_up, qds_pred_dw)

grouping = [repelem("Proliferation", 3), repelem("Quiescence", 3)]

kovatcheva_17_t = table([qds_pred_up; qds_pred_dw], grouping', 'VariableNames', {'QDS', 'grouping'})

