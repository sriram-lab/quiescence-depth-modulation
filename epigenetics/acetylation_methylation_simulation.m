%% Histone Epigenetics Simulation
% Simulates histone acetylation and methylation as a function of quiescence 
% depth by modeling substrate metabolite flux. 

%% Histone acetylation

% acetylation analysis
clear;

% loads Rat-GEM, optimized regression model, optimized params
load QDS_model_3.mat

% Rat-GEM acetylation model
model.lb(8486) = -100; % acetylated proteins 
model.ub(7137) = 1; % limited diffusion into nucleus. accoa
model.ub(7149) = 1; % limited diffusion into nucleus. coa
model.lb(3509) = 0; % choline based production
model.ub(3509) = 0; % choline based production
model.lb(11245) = 0; % citrate lyase pseudo rxn
model.ub(11245) = 0; % citrate lyase pseudo rxn

% additions: 
model.lb(3535) = 0; % choline based production - nucleus
model.ub(3535) = 0; % choline based production - nucleus

model = addReaction(model,'acl_n', 'MAM01371n + MAM01587n + MAM01597n ->	MAM01261n + MAM01285n + MAM02633n + MAM02751n');

model.grRules(end) = model.grRules(1346); % same as acly rxn
model.rxnGeneMat(end,:) = model.rxnGeneMat(1346,:);
model.rules(end) = model.rules(1346); % same as acly rxn

model.rxnNames(end) = {'ATP-Citrate lyase, nuclear'};

model = addReaction(model,'oaa_nt',  'MAM02633n -> MAM02633c');
model.lb(end) = -100;
model = addReaction(model,'cit_nt',  'MAM01587n -> MAM01587c');
model.lb(end) = -100;

model = addReaction(model,'KACYLn',' MAM01261n + MAM02708n  -> MAM03401n + MAM01597n') ;
model = addReaction(model,'peplyexn',' MAM02708n  -> MAM02708e');
model.lb(end) = -100;
model = addReaction(model,'EX_KAC',    'MAM03401n  -> ');
model = addReaction(model,'EX_KACcoa',    'MAM03401n  + MAM01597n -> ');

% addition: add peptidyl-l-lysine exchange reaction
model = addReaction(model,'pep_ex',' MAM02708e  -> ');
model.lb(end) = -100;

rxnpos  = [find(ismember(model.rxns,'EX_KAC'));];

% sets acetylation reaction objective
model.c(rxnpos) = 0.97E-3;

%% Generates fluxes across quiescence depths

% Loads quiescence deepening transcriptomics
fujimaki_19 = readtable('fujimaki_19.txt');

% isolates genes in model
fujimaki_19.Gene = upper(fujimaki_19.Gene);
fujimaki_19 = fujimaki_19(ismember(fujimaki_19.Gene, model.genes), :);

flux_q = [];

kappa = 372.7594;
rho = 372.7594;
epsilon = 0.0036;

thresh_up = 1.008695652;
thresh_down = 1 / thresh_up;


up_q = {fujimaki_19.Gene(fujimaki_19.Day3_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day3_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day3_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day4_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day4_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day4_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day6_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day6_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day6_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day8_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day8_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day8_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day10_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day10_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day10_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day12_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day12_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day12_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day14_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day14_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day14_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day16_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day16_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day16_3 ./ fujimaki_19.Day2_3 > thresh_up)
    };

dw_q = {fujimaki_19.Gene(fujimaki_19.Day3_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day3_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day3_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day4_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day4_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day4_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day6_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day6_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day6_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day8_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day8_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day8_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day10_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day10_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day10_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day12_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day12_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day12_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day14_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day14_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day14_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day16_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day16_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day16_3 ./ fujimaki_19.Day2_3 < thresh_down)
    };

%% flux generation

parfor i = 1:length(up_q)

    flux_q(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up_q{i}, dw_q{i}, kappa, rho, epsilon, 0, [], 1);
    
end

%% Acetylation flux analysis 

% takes average of triplicates for each QDS
acetylation_flux_avg = [mean(flux_q(1:3, rxnpos)), mean(flux_q(4:6, rxnpos)), ...
                 mean(flux_q(7:9, rxnpos)), mean(flux_q(10:12, rxnpos)), ...
                 mean(flux_q(13:15, rxnpos)), mean(flux_q(16:18, rxnpos)), ...
                 mean(flux_q(19:21, rxnpos)), mean(flux_q(22:24, rxnpos))];
biomass_flux_avg = [mean(flux_q(1:3, 8126)), mean(flux_q(4:6, 8126)), ...
                 mean(flux_q(7:9, 8126)), mean(flux_q(10:12, 8126)), ...
                 mean(flux_q(13:15, 8126)), mean(flux_q(16:18, 8126)), ...
                 mean(flux_q(19:21, 8126)), mean(flux_q(22:24, 8126))];

qds = [3 3 3 4 4 4 6 6 6 8 8 8 10 10 10 12 12 12 14 14 14 16 16 16]';
qds_avg = [3 4 6 8 10 12 14 16];
[r, p] = corr(acetylation_flux_avg', qds_avg')
 

clf('reset')
hold on
scatter(qds_avg, acetylation_flux_avg, 50, 'filled')
scatter(qds_avg, biomass_flux_avg, 50, 'filled')
legend(["acetylation flux", "biomass flux"])
xlabel('Quiescence depth (days since serum starvation began)')
ylabel('Mean acetylation flux')
title('Acetylation flux vs QDS (objective: 0.97E-3)')


% export to ggplot
qds_acetyl_avg_t = table(qds_avg', acetylation_flux_avg', 'VariableNames',{'qds', 'acetyl_flux_avg'});

%% Subsystems contributing to acetylation trend 

% find total subsystems in all reactions
subsystems_c_m = model.subSystems;
subsystems_m = string();

% converts cell array w/ 2 x 1 cells into string array
for i = 1:length(subsystems_c_m)

    subsystems_1_m = (subsystems_c_m{i});
    subsystems_m(i) = string(subsystems_1_m(1, 1));

end

% total unique subsystems
subsystems_unique = unique(subsystems_m);

pearson_r_subsystem_ko = zeros(length(subsystems_unique), 1);
pearson_pval_subsystem_ko = zeros(length(subsystems_unique), 1);
qds_avg = [3 4 6 8 10 12 14 16];

parfor i = 1:length(subsystems_unique)

    model_ko = model;

    % sets lower and upper bound of each reaction in subsystem to zero
    model_ko.lb(strcmp(subsystems_unique(i), subsystems_m)) = 0;
    model_ko.ub(strcmp(subsystems_unique(i), subsystems_m)) = 0;
    
    flux_q = [];

    % calculates fluxes
    for k = 1:24

        flux_q(k, :) = constrain_flux_regulation_inconsistencies(model_ko, ...
                                                      up_q{k}, dw_q{k}, kappa, rho, epsilon, 0, [], 1);
    
    end

    % average acetylation flux
    acetylation_flux_avg = [mean(flux_q(1:3, rxnpos)), mean(flux_q(4:6, rxnpos)), ...
                 mean(flux_q(7:9, rxnpos)), mean(flux_q(10:12, rxnpos)), ...
                 mean(flux_q(13:15, rxnpos)), mean(flux_q(16:18, rxnpos)), ...
                 mean(flux_q(19:21, rxnpos)), mean(flux_q(22:24, rxnpos))];

    % is there still a trend?
    [pearson_r_subsystem_ko(i), pearson_pval_subsystem_ko(i)] = corr(acetylation_flux_avg', qds_avg');

end

%% results
subsystem_ko_acetylation_t = table(subsystems_unique', pearson_r_subsystem_ko(:, 1), pearson_pval_subsystem_ko(:, 1), 'VariableNames', {'subsystem', 'pearson_r', 'pval'});
subsystem_ko_acetylation_t = sortrows(subsystem_ko_acetylation_t, 2, 'descend')



%% Methylation analysis

% loads Rat-GEM, optimized regression model, optimized params
clear;
load QDS_model_3.mat

% adds methylation reaction to GEM
model = addReaction(model,'EX_HistMET', 'MAM02129n -> '); % max nuclear histone methylation (histone-N6-methyl-L-lysine)
model.c(end) = 0.90E-3; 

% ensure methionine uptake levels positive to not limit methylation fluxes
[ix pos]  = ismember({'MAR09042'}, model.rxns);

model.lb(pos) = -0.5; % it has to be non limiting

%% Generates fluxes across quiescence depths

% Loads quiescence deepening transcriptomics
fujimaki_19 = readtable('fujimaki_19.txt');

% isolates genes in model
fujimaki_19.Gene = upper(fujimaki_19.Gene);
fujimaki_19 = fujimaki_19(ismember(fujimaki_19.Gene, model.genes), :);

flux_q = [];

kappa = 372.7594;
rho = 372.7594;
epsilon = 0.0036;

thresh_up = 1.008695652;
thresh_down = 1 / thresh_up;


up_q = {fujimaki_19.Gene(fujimaki_19.Day3_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day3_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day3_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day4_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day4_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day4_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day6_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day6_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day6_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day8_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day8_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day8_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day10_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day10_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day10_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day12_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day12_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day12_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day14_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day14_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day14_3 ./ fujimaki_19.Day2_3 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day16_1 ./ fujimaki_19.Day2_1 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day16_2 ./ fujimaki_19.Day2_2 > thresh_up), ...
    fujimaki_19.Gene(fujimaki_19.Day16_3 ./ fujimaki_19.Day2_3 > thresh_up)
    };

dw_q = {fujimaki_19.Gene(fujimaki_19.Day3_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day3_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day3_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day4_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day4_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day4_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day6_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day6_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day6_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day8_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day8_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day8_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day10_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day10_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day10_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day12_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day12_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day12_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day14_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day14_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day14_3 ./ fujimaki_19.Day2_3 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day16_1 ./ fujimaki_19.Day2_1 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day16_2 ./ fujimaki_19.Day2_2 < thresh_down), ...
    fujimaki_19.Gene(fujimaki_19.Day16_3 ./ fujimaki_19.Day2_3 < thresh_down)
    };


%% generates fluxes 

parfor i = 1:length(up_q)

    i;
    flux_q(i, :) = constrain_flux_regulation_inconsistencies(model, ...
                                                  up_q{i}, dw_q{i}, kappa, rho, epsilon, 0, [], 1);
    
end

%% Methylation flux analysis 

% takes average of triplicates for each QDS
meth_flux_avg = [mean(flux_q(1:3, 13029)), mean(flux_q(4:6, 13029)), ...
                 mean(flux_q(7:9, 13029)), mean(flux_q(10:12, 13029)), ...
                 mean(flux_q(13:15, 13029)), mean(flux_q(16:18, 13029)), ...
                 mean(flux_q(19:21, 13029)), mean(flux_q(22:24, 13029))];
biomass_flux_avg = [mean(flux_q(1:3, 8126)), mean(flux_q(4:6, 8126)), ...
                 mean(flux_q(7:9, 8126)), mean(flux_q(10:12, 8126)), ...
                 mean(flux_q(13:15, 8126)), mean(flux_q(16:18, 8126)), ...
                 mean(flux_q(19:21, 8126)), mean(flux_q(22:24, 8126))];

qds = [3 3 3 4 4 4 6 6 6 8 8 8 10 10 10 12 12 12 14 14 14 16 16 16]';
qds_avg = [3 4 6 8 10 12 14 16];
[r, p] = corr(meth_flux_avg', qds_avg')  

clf('reset')
hold on
scatter(qds_avg, meth_flux_avg, 50, 'filled')
scatter(qds_avg, biomass_flux_avg, 50, 'filled')
xlabel('Quiescence depth (days since serum starvation began)')
ylabel('Mean nuclear methylation flux')
title('Methylation flux vs QDS (objective: 0.90E-3)')
legend(["methylation flux", "biomass flux"])

qds_meth_avg_t = table(qds_avg', meth_flux_avg', 'VariableNames',{'qds', 'meth_flux_avg'})


%% Subsystem analysis 

% find total subsystems in all reactions
subsystems_c_m = model.subSystems;
subsystems_m = string()

% converts cell array w/ 2 x 1 cells into string array
for i = 1:length(subsystems_c_m)

    subsystems_1_m = (subsystems_c_m{i});
    subsystems_m(i) = string(subsystems_1_m(1, 1));

end

% total unique subsystems
subsystems_unique = unique(subsystems_m);

pearson_r_subsystem_ko = [];
pearson_pval_subsystem_ko = [];

pearson_r_subsystem_ko = zeros(length(subsystems_unique), 1);
pearson_pval_subsystem_ko = zeros(length(subsystems_unique), 1);
qds_avg = [3 4 6 8 10 12 14 16];

parfor i = 1:length(subsystems_unique)

   
    model_ko = model;

    % sets lower and upper bound of each reaction in subsystem to zero 
    model_ko.lb(strcmp(subsystems_unique(i), subsystems_m)) = 0.0;
    model_ko.ub(strcmp(subsystems_unique(i), subsystems_m)) = 0.0;

    flux_q = [];

    % calculates fluxes
    for k = 1:24

        flux_q(k, :) = constrain_flux_regulation_inconsistencies(model_ko, ...
                                                      up_q{k}, dw_q{k}, kappa, rho, epsilon, 0, [], 1);
    
    end

    % average acetylation flux
    meth_flux_avg = [mean(flux_q(1:3, 13029)), mean(flux_q(4:6, 13029)), ...
                 mean(flux_q(7:9, 13029)), mean(flux_q(10:12, 13029)), ...
                 mean(flux_q(13:15, 13029)), mean(flux_q(16:18, 13029)), ...
                 mean(flux_q(19:21, 13029)), mean(flux_q(22:24, 13029))];

    % is there still a trend?
    [pearson_r_subsystem_ko(i), pearson_pval_subsystem_ko(i)] = corr(meth_flux_avg', qds_avg');

end

%% results
subsystem_ko_methylation_t = table(subsystems_unique', pearson_r_subsystem_ko(:, 1), pearson_pval_subsystem_ko(:, 1), 'VariableNames', {'subsystem', 'pearson_r', 'pval'});

subsystem_ko_methylation_t = sortrows(subsystem_ko_methylation_t, 2)