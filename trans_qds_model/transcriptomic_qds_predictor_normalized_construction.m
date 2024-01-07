%% Builds a QDS predictor based on multiple transcriptomics datasets

% Builds and optimizes a multi-cell-type quiescence depth predictor using
% gene expression levels. Predicts qds value normalized for each dataset 
% to a value between [0, 1]. 

clear

%% loads transcriptomics + qds csvs

% tpm + quantile normalized
transcriptomics = readtable("aggregated_trans_quantile_normalized.csv")

genes = transcriptomics.Var1

transcriptomics_mat = table2array(transcriptomics(:, 2:end));

qds_normalized = readmatrix("aggregated_qds_1_normalized.csv")

%% Optimizing alpha (ratio of ridge to lasso) 

% wide sweep
alpha_array = [];
alpha_array_no_repeats = logspace(-4, 0, 500);

alpha_array = repelem(alpha_array_no_repeats, 3);

mse_cv = [];
num_coeff = [];
min_lambda = [];

X_train = []; Y_train = [];
X_train = transcriptomics_mat';
Y_train = qds_normalized;

% z-score normalize x_train
X_train = zscore(X_train, 1);

parpool(7)

parfor i = 1:length(alpha_array)

    i
    B = []; FitInfo = [];

    % estimates MSE through 10-fold cross-validation
    [B, FitInfo] = lasso(X_train, Y_train, 'Alpha', alpha_array(i), 'CV', 10, 'Standardize', true);

    min_lambda(i) = FitInfo.LambdaMinMSE

    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    coef = B(:, idxLambdaMinMSE);
     
    % finds crossfold mse
    mse_cv(i) = FitInfo.MSE(idxLambdaMinMSE)

    num_coeff(i) = sum(coef ~= 0)

end

%% takes average of 3 replicates of mse
mse_cv_avg = [];
num_coeff_avg = [];
min_lambda_avg = [];

k = 1

for i = 1:3:1500

    mse_cv_avg(k) = mean(mse_cv(i:(i+2)));
    num_coeff_avg(k) = mean(num_coeff(i:(i+2)));
    min_lambda_avg(k) = mean(min_lambda(i:(i+2)));

    k = k + 1;

end

plotyy(log10(alpha_array_no_repeats), (mse_cv_avg), log10(alpha_array_no_repeats), num_coeff_avg)
legend(["CV MSE", "Num Coeff"])
xlabel("Alpha")

[M, I] = min(mse_cv_avg)

alpha_array_no_repeats(I)
num_coeff_avg(I)

%% narrow sweep

alpha_array = [];

% range: optimal from wide sweep - 0.003 to optimal from wide sweep + 0.003
alpha_array_no_repeats = linspace(0.0009, 0.0069, 500);

alpha_array = repelem(alpha_array_no_repeats, 3);

mse_cv = [];
num_coeff = [];

X_train = []; Y_train = [];
X_train = transcriptomics_mat';
Y_train = qds_normalized;

% z-score normalize x_train
X_train = zscore(X_train, 1);

parpool(7)

parfor i = 1:length(alpha_array)

    i
    B = []; FitInfo = [];

    [B, FitInfo] = lasso(X_train, Y_train, 'Alpha', alpha_array(i), 'CV', 10, 'Standardize', true);

    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    coef = B(:, idxLambdaMinMSE);
    
    % finds crossfold mse
    mse_cv(i) = FitInfo.MSE(idxLambdaMinMSE)
    num_coeff(i) = sum(coef ~= 0)

end

%% takes average of 3 replicates of mse
mse_cv_avg = [];
num_coeff_avg = [];

k = 1

for i = 1:3:1500

    mse_cv_avg(k) = mean(mse_cv(i:(i+2)));
    num_coeff_avg(k) = mean(num_coeff(i:(i+2)));

    k = k + 1

end

plotyy(log10(alpha_array_no_repeats), (mse_cv_avg), log10(alpha_array_no_repeats), num_coeff_avg)
legend(["CV MSE", "Num Coeff"])
xlabel("Alpha")

[M, I] = min(mse_cv_avg)

alpha_array_no_repeats(I)
num_coeff_avg(I)

%% Building regression model with optimal parameters

% all
X_train = transcriptomics_mat';
Y_train = qds_normalized;

% z-score normalize x_train
X_train = zscore(X_train, 1);

alpha_optimal = 0.002378958;


[B, FitInfo] = lasso(X_train, Y_train, 'Alpha', alpha_optimal, 'CV', 10, 'Standardize', true);

idxLambdaMinMSE = FitInfo.IndexMinMSE
coef = B(:, idxLambdaMinMSE);
coef0 = FitInfo.Intercept(idxLambdaMinMSE);

yhat = X_train * coef + coef0;

lambda_min = FitInfo.LambdaMinMSE

sum(coef ~= 0)

mse = FitInfo.MSE(idxLambdaMinMSE)


%% Saving of optimized linear regression model 

% creates coefficients + genes table
coef_table = table(transcriptomics.Var1, coef)

coef_table.Properties.VariableNames = {'gene', 'coef'}

% for later use in ggplot
genes_ggplot = transcriptomics.Var1


%% Random index 10-fold cross validation 

X_train = []; Y_train = [];
X_train = transcriptomics_mat';
Y_train = qds_normalized;

cv_idx = crossvalind('kfold', Y_train, 10)

X_test_predictions = zeros([length(Y_train) 1]);

for i = 1:10

    i
    idx_train = cv_idx ~= i;
    idx_test = cv_idx == i;

    X_train_subset = X_train(idx_train, :);
    Y_train_subset = Y_train(idx_train);

    X_test_subset = X_train(idx_test, :);
    Y_test_subset = Y_train(idx_test);

    % z-score normalize X_train and X_test
    X_train_subset = zscore(X_train_subset, 1);
    X_test_subset = zscore(X_test_subset, 1);


    % optimize alpha: 
    alpha_array = linspace(1E-4, 1, 100)

    mse_cv = [];
    num_coeff = [];

    parfor j = 1:length(alpha_array)
            
        B = []; FitInfo = [];
    
        [B, FitInfo] = lasso(X_train, Y_train, 'Alpha', alpha_array(j), 'CV', 10, 'Standardize', true);
    
        idxLambdaMinMSE = FitInfo.IndexMinMSE;
        
        % finds crossfold mse
        mse_cv(j) = FitInfo.MSE(idxLambdaMinMSE)

    end


    [M, I] = min(mse_cv)
    
    min_alpha = alpha_array(I)

    % builds model with optimal parameters
    [B, FitInfo] = lasso(X_train_subset, Y_train_subset, 'Alpha', min_alpha, 'CV', 10, 'Standardize', true);
    
    idxLambdaMinMSE = FitInfo.IndexMinMSE
    coef = B(:, idxLambdaMinMSE);
    coef0 = FitInfo.Intercept(idxLambdaMinMSE);
        
    lambda_min = FitInfo.LambdaMinMSE

    % applies model to data it's never seen before
    X_test_predictions(idx_test) = X_test_subset * coef + coef0;

end

%% Evaluate model accuracy

scatter(Y_train, X_test_predictions)

[r, p] = corr(Y_train, X_test_predictions)
