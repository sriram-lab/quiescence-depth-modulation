# Normalizes, consolidates, and omits outliers for transcriptomics datasets used to build
# multi-cell-type QDS predictor


### install packages
install.packages("devtools")
devtools::install_version("dbplyr", version = "2.3.4")

### load packages
library(tidyverse)
library(biomaRt)
library(datawizard)
library(gridExtra)
library(preprocessCore)
library(testit)
library(factoextra)
library(gtools)
library(devtools)
library(dbplyr)

# define path to input files
my_directory <- ""

### sturm_23
### already in TPM, filter to just include contact inhibition

sturm_23 <- read.csv(file.path(my_directory, "datasets", "sturm_23", "sturm_23_GSE179848_txi_gene_count_lengthScaledTPM.csv"))

sturm_23 <- as.data.frame(sturm_23)

# convert gene name to "gene"
colnames(sturm_23)[1] <- "gene"

# sort by columns by numeric value so it lines up with GEO
sturm_23 <- sturm_23[, mixedsort(colnames(sturm_23))]

# contact inhibition idx: sample 217 to 257, column 218 to 258
# 41 samples

sturm_23 <- sturm_23[, c(1, 218:258)]

sturm_23_tpm <- sturm_23

sturm_23_qds <- c(rep(seq(0, 120, 20), 5), 0, 20, 60, 80, 100, 120)

# omit day 0 (since most likely proliferating here)
idx_omit_sturm_23 <- which(sturm_23_qds == 0)

sturm_23_qds <- sturm_23_qds[-idx_omit_sturm_23]

# omit same indices but add 1 to account for gene column
sturm_23_tpm <- sturm_23_tpm[, -(idx_omit_sturm_23 + 1)]


### fujimaki_19
# convert from FPKM to TPM

fujimaki_19 <- read.table(file.path(my_directory, "datasets", "fujimaki_19_fpkm.txt"), header = T) %>% as.data.frame()

# convert gene name to "gene"
colnames(fujimaki_19)[1] <- "gene"

fujimaki_19[, 1] <- toupper(fujimaki_19[, 1])

fujimaki_19_fpkm = fujimaki_19[, 2:ncol(fujimaki_19)]

fujimaki_19_tpm <- as.data.frame(apply(fujimaki_19_fpkm, 2, fpkm_to_tpm))

# omit controls
fujimaki_19_tpm <- fujimaki_19_tpm[, 4:ncol(fujimaki_19_tpm)]

fujimaki_19_tpm$gene <- fujimaki_19$gene

fujimaki_19_qds <- c(rep(2, 3), rep(3, 3), rep(4, 3), rep(6, 3), 
                     rep(8, 3), rep(10, 3), rep(12, 3), rep(14, 3), rep(16, 3))


### rajapakse_20

# convert from counts to TPM

rajapakse_20 = read.table(file.path(my_directory, "datasets", "rajapakse_20_GSE129964_serumCountsTable.tsv"), header = T, row.names = 1)

# add gene column
rajapakse_20$gene <- rownames(rajapakse_20)

# omit proliferating controls
rajapakse_20 <- rajapakse_20[, 4:22]

# filter to only genes in fujimaki_19
rajapakse_20 <- rajapakse_20[rajapakse_20$gene %in% fujimaki_19$gene, ]

# filter out genes with low counts: all gene for all conditions must have count ≥ 10
count_filter <- rowSums(rajapakse_20[, 1:18] >= 10) >= (ncol(rajapakse_20) - 1)

rajapakse_20 <- rajapakse_20[count_filter, ]


# get gene lengths
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

transcript_length_info <- getBM(attributes = c("ensembl_gene_id", "transcript_length", "external_gene_name"), values = rajapakse_20$gene, mart = ensembl)

transcript_length_mean <- transcript_length_info %>% group_by(external_gene_name) %>% summarise(transcript_length_avg = mean(transcript_length))

colnames(transcript_length_mean)[1] <- "gene"

# merge with rajapakse_20
rajapakse_20_with_transcript_length <- left_join(rajapakse_20, transcript_length_mean, by = "gene")

# omit nas
rajapakse_20_with_transcript_length <- rajapakse_20_with_transcript_length %>% na.omit()


rajapakse_20_counts <- rajapakse_20_with_transcript_length[, 1:18]
rajapakse_20_lengths <- rajapakse_20_with_transcript_length[, 20] %>% as.numeric()

rajapakse_20_tpm <- matrix(0, nrow = nrow(rajapakse_20_counts), ncol = ncol(rajapakse_20_counts))


for (i in 1:ncol(rajapakse_20_counts)) {

  rajapakse_20_tpm[, i] <- count_to_tpm(rajapakse_20_counts[, i], rajapakse_20_lengths)
}

rajapakse_20_tpm <- rajapakse_20_tpm %>% as.data.frame()

colnames(rajapakse_20_tpm) <- colnames(rajapakse_20_counts)

rajapakse_20_tpm$gene <- toupper(rajapakse_20_with_transcript_length$gene)


rajapakse_20_qds <- c(c(rep(1, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3), rep(9, 3)))



# aggregate all TPMs into one dataframe, merge by gene symbol ("gene")

df_list <- list(fujimaki_19_tpm, rajapakse_20_tpm, sturm_23_tpm)

# un-normalized qds
aggregated_qds <- c(fujimaki_19_qds, rajapakse_20_qds, sturm_23_qds)

# normalized qds (each qds normalized to [0 1])
aggregated_qds_normalized <- c(fujimaki_19_qds / 16, rajapakse_20_qds / 9, sturm_23_qds / 120)


# distribution of qds values
hist(aggregated_qds)

# merge all data frames in list
aggregated_trans <- df_list %>% purrr::reduce(inner_join, by='gene')

# omit duplicate gene rows
aggregated_trans <- aggregated_trans[!duplicated(aggregated_trans$gene), ]

# move gene column to row names
rownames(aggregated_trans) <- aggregated_trans$gene

# remove gene column
aggregated_trans <- aggregated_trans[, -(which(colnames(aggregated_trans) == "gene"))]

aggregated_trans <- mutate_all(aggregated_trans, as.numeric)

# save column labels and gene names
aggregated_labels <- colnames(aggregated_trans)
aggregated_genes <- rownames(aggregated_trans)

# run quantile normalization to normalize sample distribution

data_source_color <- c(rep(2, 27), rep(3, 18), rep(1, 35))

png(file.path(my_directory, "trans_qds_predictor", "results", "plots", "transcript_distribution_before.png"), res = 350, width = 5, height = 5, units = "in")

# plot distribution of expression values before normalization
plot(density(aggregated_trans[, 1]), ylim = c(0, 0.04), xlim = c(0, 2000), col = data_source_color[1], xlab = "Transcript level", main = "Transcript Distribution:\nBefore Quantile Normalization")

for (i in 2:ncol(aggregated_trans)) {
  
  lines(density(aggregated_trans[, i]), add = T, col = data_source_color[i])
}

legend(800, 0.035, legend=c("REF: GSE124109", "RPE: GSE129964", "HFB: GSE179848"),
       col=c(2, 3, 1), lty=1:1, cex=1,
       title="", text.font=4, bg='white')

dev.off()

# quantile normalize
aggregated_trans <- as.data.frame(normalize.quantiles(as.matrix(aggregated_trans)))

# add row and column names back
colnames(aggregated_trans) <- aggregated_labels
rownames(aggregated_trans) <- aggregated_genes

png(file.path(my_directory, "trans_qds_predictor", "results", "plots", "transcript_distribution_after_quantile_norm.png"), res = 350, width = 5, height = 5, units = "in")

# check distribution after, seems pretty good
plot(density(aggregated_trans[, 1]), ylim = c(0, 0.04), xlim = c(0, 2000), col = data_source_color[1], xlab = "Transcript level", main = "Transcript Distribution:\nAfter Quantile Normalization")

for (i in 2:ncol(aggregated_trans)) {
  
  lines(density(aggregated_trans[, i]), add = T, col = data_source_color[i])
}

legend(800, 0.035, legend=c("REF: GSE124109", "RPE: GSE129964", "HFB: GSE179848"),
       col=c(2, 3, 1), lty=1:1, cex=1,
       title="", text.font=4, bg='white')

dev.off()

# run PCA + hierarchical clustering + more QC 

res.pca <- prcomp(t(aggregated_trans), scale = F)

fviz_eig(res.pca)


data_source_color <- c(rep("REF: GSE124109", 27), rep("RPE: GSE129964", 18), rep("HFB: GSE179848", 35))

pca_before <- fviz_pca_ind(res.pca,
             habillage = data_source_color,
             repel = TRUE,    
             labelsize = 1,
             addEllipses = T,
             pointsize = 1,
             label = "none",
             palette = c("black", "red", "green3")
) + labs(title = "PCA: Before Filtering")

print(pca_before)

ggsave(file.path(my_directory, "trans_qds_predictor", "results", "plots", "pca_before_outlier_removal.png"), plot = pca_before, dpi = 350, width = 4, height = 4.5)

# outliers: s.day4.repeat2
outlier_idx <- c(35)

aggregated_trans <- aggregated_trans[, -outlier_idx]
aggregated_qds <- aggregated_qds[-outlier_idx]
aggregated_qds_normalized <- aggregated_qds_normalized[-outlier_idx]

# run pca after outliers removed to double check they're not there
res.pca <- prcomp(t(aggregated_trans), scale = F)

fviz_eig(res.pca)

pca_after <- fviz_pca_ind(res.pca,
                           habillage = data_source_color[-outlier_idx],
                           repel = TRUE,     
                           labelsize = 1,
                           addEllipses = T,
                           pointsize = 1,
                           label = "none",
                           palette = c("black", "red", "green3")
) + labs(title = "PCA: After Filtering")

print(pca_after)

ggsave(file.path(my_directory, "trans_qds_predictor", "results", "plots", "pca_after_outlier_removal.png"), plot = pca_after, dpi = 350, width = 4, height = 4.5)

# cluster by condition / cell type, as expected

# z-score normalize
aggregated_trans_z_score <- t(scale(t(aggregated_trans)))

write.csv(aggregated_trans_z_score, file.path(my_directory, "aggregated_trans_z_score_quantile_normalized.csv"), row.names = T)

write.csv(aggregated_qds_normalized, file.path(my_directory, "aggregated_qds_1_normalized.csv"), row.names = F)


## plotting distribution of QDS and highlighting data type
qds_df <- data.frame(
  qds_norm = aggregated_qds_normalized,
  qds = aggregated_qds,
  data_source = data_source_color[-outlier_idx]
)

qds_density_unnormalized <- ggplot(qds_df, aes(x = qds, fill = data_source, color = data_source)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(name = "", values = c("black", "red", "green3")) +
  scale_color_manual(name = "", values = c("black", "red", "green3")) +
  labs(x = "QDS", y = "Density", title = "QDS Distribution: Unnormalized") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15))
qds_density_unnormalized

ggsave(file.path(my_directory, "trans_qds_predictor", "results", "plots", "qds_distribution_unnormalized.png"), plot = qds_density_unnormalized, dpi = 350, width = 4, height = 4.5)

qds_density_normalized <- ggplot(qds_df, aes(x = qds_norm, fill = data_source, color = data_source)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(name = "", values = c("black", "red", "green3")) +
  scale_color_manual(name = "", values = c("black", "red", "green3")) +
  labs(x = "QDS (normalized)", y = "Density", title = "QDS Distribution: Normalized") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5))
qds_density_normalized

ggsave(file.path(my_directory, "trans_qds_predictor", "results", "plots", "qds_distribution_normalized.png"), plot = qds_density_normalized, dpi = 350, width = 4, height = 4.5)



##### FUNCTIONS ######

# converts from counts to tpm
count_to_tpm <- function(counts, length)
{
  numerator <- counts / length
  denominator <- sum(numerator)

  tpm <- (numerator / denominator) * 1e6
}

# converts from fpkm to tpm
fpkm_to_tpm <- function(fpkm)
{
  tpm <- (fpkm / sum(fpkm)) * 1e6
}

