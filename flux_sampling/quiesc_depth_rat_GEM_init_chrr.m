%% Predict metabolic changes underlying quiescence deepening with INIT + CHRR

% For each stage of quiescence depth (from 2 to 16 days of serum
% starvation), transcriptomics-constrained model is created with INIT, then
% fluxes are generated through sampling with CHRR

%% Initializing COBRA 

clear;

% initializes cobra toolbox + gurobi solver
initCobraToolbox(false);
changeCobraSolver('gurobi');

%% Loading quiescence deepening transcriptomics
fujimaki_19 = readtable('fujimaki_19.txt');


%% Loading Rat-GEM
load Rat-GEM.mat;
model = ratGEM;

% generates rules
model.b = model.b(:, 1);
model = generateRules(model);

% converts gene symbols to uppercase
model.genes = upper(model.genes);

% sets objective function to 0
model.c(model.c == 1) = 0;


%% Further processing of transcriptomics

% isolates genes in Rat-GEM metabolic model
fujimaki_19.Gene = upper(fujimaki_19.Gene);
fujimaki_19 = fujimaki_19(ismember(fujimaki_19.Gene, model.genes), :);

% isolates transcript levels
fujimaki_19_mat = table2array(fujimaki_19(:, 2:end));

% quantile normalizes
norm_data = quantilenorm(fujimaki_19_mat);

% log2 normalizes
norm_data = log2(norm_data + 1);


%% compute fluxes for each stage of quiescence

% sets indices of each quiescence depth, from 2 to 16 (skipping first 3,
% proliferation)
idx_qds = {4:6 7:9 10:12 13:15 16:18 19:21 22:24 25:27 28:30};

% creates empty cell array for constrained models
quiescence_models = {};

% runs INIT on each quiescence depth to generate series of constrained
% models
for i = 1:length(idx_qds)

    disp(strcat("Processing ", string(i), " of ", string(length(idx_qds)), " samples"))

    norm_data_subset = norm_data(:, idx_qds{i});
    
    alpha = 1.1; beta = 0.6;

    % runs INIT to obtain metabolic model constrained by transcriptomics
    quiescence_models{i} = run_init2(model, alpha, beta, norm_data_subset, fujimaki_19.Gene, 32);

end

flux_avg_array = {};
flux_all_array = {};

% creates parallel pool
pool = parpool('Processes', 9)

environment = getEnvironment();

% runs CHRR on each quiescence depth 
parfor i = 1:length(quiescence_models)

    disp(i)
    num_samples = 1000

    restoreEnvironment(environment)

    % runs CHRR
    [flux_avg_array{i}, flux_all_array{i}] = run_chrr2(model, quiescence_models{i}, num_samples);

end

delete(pool)
