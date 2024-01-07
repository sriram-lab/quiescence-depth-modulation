%% Embryogenesis wang_21 corrobration

% Looks at expression levels of proliferation-related genes to corroborate
% fall + rise in QDS

% loads imputed embryogenesis data 
wang_21 = readtable("wang_21_embryogenesis_imputed_filtered_go_prolif_genes.csv")

wang_21.gene = upper(wang_21.gene);

counts = wang_21{:, 2:end};

% Normalization: scaling each column by size factor
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

normCounts = bsxfun(@rdivide,counts,sizeFactors);
normCounts(1:10,:);

wang_21{:, 2:end} = normCounts;

% omit missing values
wang_21 = rmmissing(wang_21);

%% loads GO gene sets

% loads GO cell growth genes
go_cell_growth_t = readtable("GOBP_CELL_GROWTH.v2023.2.Mm.csv")

go_cell_growth_genes = upper(split(go_cell_growth_t{17, 2}, ','))

% loads stem cell prolif genes
go_sc_prolif_t = readtable("GOBP_STEM_CELL_PROLIFERATION.v2023.2.Hs.csv")

go_sc_prolif_genes = upper(split(go_sc_prolif_t{17, 2}, ','))


%% GO cell growth gene set

wang_21_cell_growth = wang_21(ismember(wang_21.gene, go_cell_growth_genes), :);

wang_21_cell_growth_mat = wang_21_cell_growth{:, 2:end};

% z-score normalize
wang_21_cell_growth_mat = normalize(wang_21_cell_growth_mat, 2);

wang_21_cell_growth_mat_avg = mean(wang_21_cell_growth_mat, 1)

grouping = strings()
grouping(2:27) = "Zygote"
grouping(28:75) = "2 Cell"
grouping(76:169) = "4 Cell"
grouping(170:217) = "Late 4 Cell"
grouping(218:281) = "8 Cell"
grouping(282:329) = "16 Cell"
grouping(330:402) = "32 Cell"

grouping = grouping(2:end)

p = anova1(wang_21_cell_growth_mat_avg, grouping)

embryo_cell_growth_t = table(wang_21_cell_growth_mat_avg', grouping', 'VariableNames', {'go_bp_cell_growth_zscore_avg', 'grouping'})


%% GO stem cell proliferation gene set

wang_21_sc_prolif = wang_21(ismember(wang_21.gene, go_sc_prolif_genes), :);

wang_21_sc_prolif_mat = wang_21_sc_prolif{:, 2:end};

% z-score normalize
wang_21_sc_prolif_mat = normalize(wang_21_sc_prolif_mat, 2);

wang_21_sc_prolif_mat_avg = mean(wang_21_sc_prolif_mat, 1)

grouping = strings()
grouping(2:27) = "Zygote"
grouping(28:75) = "2 Cell"
grouping(76:169) = "4 Cell"
grouping(170:217) = "Late 4 Cell"
grouping(218:281) = "8 Cell"
grouping(282:329) = "16 Cell"
grouping(330:402) = "32 Cell"

grouping = grouping(2:end)

p = anova1(wang_21_sc_prolif_mat_avg, grouping)

embryo_sc_prolif_t = table(wang_21_sc_prolif_mat_avg', grouping', 'VariableNames', {'zscore_avg', 'grouping'})
