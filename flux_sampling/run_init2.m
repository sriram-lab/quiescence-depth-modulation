% Given initial metabolic model and omics input, returns constrained model
% based on INIT

% INPUTS:
% initial_model: unconstrained metabolic model

% alpha: upregulated genes must be greater than mean + alpha * standard
% deviation, downregulated genes must be less than mean - alpha * standard
% deviation

% beta: at least beta percent of genes must meet the criteria described
% just above

% norm_trans: quantile-normalized, log2 normalized omics data

% gene_list: list of genes corresponding to rows in norm_trans and aligning
% with initial_model gene format

% num_threads: number of threads for parallel processing

% OUTPUTS:
% constrained_model: metabolic model constrained by omics gene levels

function [constrained_model] = run_init2(initial_model, alpha, beta, norm_trans, gene_list, num_threads)


    % ensures gene list and number of rows in gene list are the same size
    assert(length(gene_list) == height(norm_trans))

    % finds highly and lowly expressed genes in transcriptomics data
    
    mean_trans = mean(norm_trans);
    stdev_trans = std(norm_trans);
    
    idx_high = sum(norm_trans > mean_trans + alpha .* stdev_trans, 2) >= beta .* width(norm_trans);
    
    disp(strcat("Number upregulated genes: ", string(sum(idx_high))))
    
    idx_low = sum(norm_trans < mean_trans - alpha .* stdev_trans, 2) >= beta .* width(norm_trans);
    
    disp(strcat("Number downregulated genes: ", string(sum(idx_low))))
    
    % converts gene expression to format for INIT input
    
    gene_expression_input = zeros([height(norm_trans) 1]);
    
    gene_expression_input(idx_high) = 1;
    gene_expression_input(idx_low) = -1;
    
    % multiply by 2 to avoid -1 input, which can cause problems
    gene_expression_input = gene_expression_input .* 2;
    
    % prepare data for INIT
    trans_data = struct();
    trans_data.gene = gene_list;
    trans_data.value = gene_expression_input;
    
    trans_rxn_expression = mapExpressionToReactions(initial_model, trans_data);
    
    % sets INIT options
    options = struct();
    options.solver = 'INIT';   
    options.weights = trans_rxn_expression;
    options.timelimit = 45;
    options.printLevel = 1;
    options.numThreads = num_threads;
 
    % runs INIT
    constrained_model = createTissueSpecificModel(initial_model, options, 1)

end
