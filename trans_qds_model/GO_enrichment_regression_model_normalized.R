# Determines GO enrichment for positive and negative coefficients in normalized 
# regression model out of all genes used to train the model (background)

### intall packages
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
install.packages("R.matlab")
BiocManager::install("biocLite.R")

### load packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(R.matlab)
library(DOSE)

### define directories
my_directory <- ""
my_directory_datasets <- ""
my_directory_datasets_tpm <- ""
plots_directory <- file.path("")
data_directory_trans_qds <- file.path("")


## read in regression model data
regression_model_normalized <- readMat(file.path(my_directory, "multi_cell_type_transcriptomics_quiescence_depth_predictor.mat"))

coefs <- regression_model_normalized$coef

genes <- unlist(regression_model_normalized$genes)

genes <- (regression_model_normalized$genes) %>% unlist() %>% as.data.frame()

## read in training data: genes here are background set
training_data <- read.csv(file.path(my_directory, "aggregated_trans_z_score_quantile_normalized.csv"))


background_gene_set <- training_data$X

## perform GO enrichment - up
target_genes_up <- genes[coefs > 0]


GO_results_up <- enrichGO(gene = target_genes_up, 
                          universe = background_gene_set,
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "SYMBOL", 
                          ont = "BP",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05) %>% as.data.frame() 

GO_results_up$fold_enrichment <- parse_ratio(GO_results_up$GeneRatio) / parse_ratio(GO_results_up$BgRatio)

GO_results_up<- GO_results_up %>% na.omit()

write.csv(GO_results_up, file.path(data_directory_trans_qds, "normalized_regression_model_pos_coefficients_go_bp_all_qvalue_0.05.csv"))

# show top 8
GO_results_up <- GO_results_up[1:8, ]


qds_go_lm_model_up_plot <- ggplot(GO_results_up, aes(x = reorder(Description, fold_enrichment), y = fold_enrichment, fill = qvalue)) +
  geom_point(pch = 21, stroke = 0, size = 2) +
  scale_fill_gradient(name = "Q-value", high = "#75ddff", low = "navy") +
  labs(x = "", y = "Fold enrichment", title = "Up in Quiescence Deepening") +
  #scale_y_continuous(breaks = c(5, 15, 25)) +
  coord_flip() +
  theme(legend.margin = margin(l = 10), 
        legend.spacing.y = unit(0.4, 'cm'),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 8),
        panel.background = element_rect(fill='grey94', colour='white'), 
        plot.title = element_text(hjust = 1))
print(qds_go_lm_model_up_plot)

ggsave(file.path(plots_directory, "normalized_regression_model_pos_coefficients_go_bp_top_8_qvalue_0.05.png"), 
       plot = qds_go_lm_model_up_plot, dpi = 400, scale = 2.9, width = 2, height = 0.9, units = "in")



## perform GO enrichment - down
target_genes_dw <- genes[coefs < 0]

GO_results_dw <- enrichGO(gene = target_genes_dw, 
                          universe = background_gene_set,
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "SYMBOL", 
                          ont = "BP",
                          pvalueCutoff = 0.65,
                          qvalueCutoff = 0.65) %>% as.data.frame() 


GO_results_dw$fold_enrichment <- parse_ratio(GO_results_dw$GeneRatio) / parse_ratio(GO_results_dw$BgRatio)

GO_results_dw <- GO_results_dw %>% na.omit()

write.csv(GO_results_dw, file.path(data_directory_trans_qds, "normalized_regression_model_neg_coefficients_go_bp_all_qvalue_0.65.csv"))

# show top 8
GO_results_dw <- GO_results_dw[1:8, ]

qds_go_lm_model_dw_plot <- ggplot(GO_results_dw, aes(x = reorder(Description, fold_enrichment), y = fold_enrichment, fill = qvalue)) +
  geom_point(pch = 21, stroke = 0, size = 2) +
  scale_fill_gradient(name = "Q-value", high = "#75ddff", low = "navy") +
  labs(x = "", y = "Fold enrichment", title = "Down in Quiescence Deepening") +
  #scale_y_continuous(breaks = c(5, 15, 25)) +
  coord_flip() +
  theme(legend.margin = margin(l = 10), 
        legend.spacing.y = unit(0.4, 'cm'),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 8),
        panel.background = element_rect(fill='grey94', colour='white'), 
        plot.title = element_text(hjust = 1))
print(qds_go_lm_model_dw_plot)

ggsave(file.path(plots_directory, "normalized_regression_model_neg_coefficients_go_bp_top8_qvalue_0.65.png"), 
       plot = qds_go_lm_model_dw_plot, dpi = 400, scale = 2.6, width = 2, height = 1, units = "in")

