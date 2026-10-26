library(TAILcaller)
library(Rsamtools)
library(dplyr)
library(stats)
library(tidyr)
library(rlang)
library(ggplot2)
library(ggrastr)
library(ComplexHeatmap)
library(ggtree)
library(circlize)
library(car)
library(dunn.test)
library(nortest)
library(rstatix)
library(rtracklayer)
library(stringr)
library(viridis)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd('C:/Users/bacs5867/Documents/PABPN1/scripts/notebooks_publi')

#############################################################
#heatmap

df <- read.table("Fig1/data_heatmap_polya.txt",header = TRUE,check.names = FALSE)
row.names(df) <- df$PolyA_lengths
df <- df[,-1]
df[] <- lapply(df, function(x) as.numeric(as.character(x)))


mat <- t(as.matrix(df))
rownames(mat)[1:5]
colnames(mat)[1:5]


order_rows <- c(
  "shCTL_nucleus_1",
  "shCTL_nucleus_2",
  "shCTL_nucleus_3",
  "shPABPN1_nucleus_1",
  "shPABPN1_nucleus_2",
  "shPABPN1_nucleus_3",
  "shCTL_cytoplasm_1",
  "shCTL_cytoplasm_2",
  "shCTL_cytoplasm_3",
  "shPABPN1_cytoplasm_1",
  "shPABPN1_cytoplasm_2",
  "shPABPN1_cytoplasm_3"
)

order_columns <- colnames(mat)

mat_df <- as.data.frame(mat) %>%
  rownames_to_column("Condition") %>%
  pivot_longer(
    -Condition,
    names_to = "Bin",
    values_to = "Value"
  )

mat_df <- mat_df %>%
  mutate(
    Condition = factor(Condition, levels = rev(order_rows)),
    Bin = factor(Bin, levels = order_columns)
  )


p <- ggplot(mat_df, aes(x = Bin, y = Condition, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#1B9E77",
    mid = "white",
    high = "#D95F02",
    midpoint = median(mat_df$Value),
    name = "Reads\nproportion (%)"
  ) +
  labs(x = "PolyA length (bp)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
p


#############################################################
# volcano plots poly(A) length

nuc_genes_polya_compare <- read.table('Fig1/nuc_genes_polya_compare.txt', sep = '\t', header = T)
nuc_genes_polya_compare$X.log10pvalue <- -log10(nuc_genes_polya_compare$FDR)

volcano <- ggplot(data = nuc_genes_polya_compare, aes(log2FoldChange, X.log10pvalue, col = effect.on.polyA.tail)) + 
  geom_point_rast(size = 2, show.legend = FALSE) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  #ggtitle('Nucleus')+
  ylab(expression(-log[10]~"(p-value)"))+
  xlab(expression(log[2]~Fold~Change~"(median poly(A) tail length)"))+
  xlim(-2.5, 2.5)+
  ylim(0, 350) +
  scale_color_manual(values = c("grey", "#984EA3", "#4DAF4A"), labels = c("", "", ""), name = '')+
  theme_bw(base_size = 20)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano


cyto_genes_polya_compare <- read.table('Fig1/cyto_genes_polya_compare.txt', sep = '\t', header = T)
cyto_genes_polya_compare$X.log10pvalue <- -log10(cyto_genes_polya_compare$FDR)

volcano <- ggplot(data = cyto_genes_polya_compare, aes(log2FoldChange, X.log10pvalue, col = effect.on.polyA.tail)) + 
  geom_point_rast(size = 2, show.legend = FALSE) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  #ggtitle('Cytoplasm')+
  ylab(expression(-log[10]~"(p-value)"))+
  xlab(expression(log[2]~Fold~Change~"(median poly(A) tail length)"))+
  xlim(-2.5, 2.5)+
  ylim(0, 350) +
  scale_color_manual(values = c("grey", "#984EA3", "#4DAF4A"), labels = c("", "", ""), name = '')+
  theme_bw(base_size = 20)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano


#############################################################
# GO analysis - cluster enrichments

genes_polya <- read.table('Fig1/genes_polya_median_clusters.txt', sep = '\t', header = T)
background_genes <- genes_polya$ensembl_gene # genes with 20 reads in each condition

genes_cluster1 <- as.list(genes_polya[genes_polya$cluster == 1,]$ensembl_gene)
genes_cluster2 <- as.list(genes_polya[genes_polya$cluster == 2,]$ensembl_gene)
genes_cluster3 <- as.list(genes_polya[genes_polya$cluster == 3,]$ensembl_gene)
genes_cluster4 <- as.list(genes_polya[genes_polya$cluster == 4,]$ensembl_gene)
genes_cluster5 <- as.list(genes_polya[genes_polya$cluster == 5,]$ensembl_gene)

#########
GO_results_cluster1 = enrichGO(gene = genes_cluster1, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'CC', pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                               universe = background_genes)
res_df_cluster1 <- as.data.frame(GO_results_cluster1)
########
GO_results_cluster2 = enrichGO(gene = genes_cluster2, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'CC', pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                               universe = background_genes)
res_df_cluster2 <- as.data.frame(GO_results_cluster2)
########
GO_results_cluster3 = enrichGO(gene = genes_cluster3, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'CC', pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                               universe = background_genes)
res_df_cluster3 <- as.data.frame(GO_results_cluster3)
########
GO_results_cluster4 = enrichGO(gene = genes_cluster4, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'CC', pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                               universe = background_genes)
res_df_cluster4 <- as.data.frame(GO_results_cluster4)
########
GO_results_cluster5 = enrichGO(gene = genes_cluster5, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'CC', pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                               universe = background_genes)
res_df_cluster5 <- as.data.frame(GO_results_cluster5)
########

# plot cluster 2 enrichments
terms_keep <- c('ribosome', 'mitochondrial inner membrane', 'mitochondrial matrix', 'nuclear pore', 'chromosomal region') # remove redundant terms
idx <- GO_results_cluster2@result$Description %in% terms_keep
GO_results_cluster2_selected <- GO_results_cluster2
GO_results_cluster2_selected@result <- GO_results_cluster2@result[idx, ]

p <- dotplot(GO_results_cluster2_selected, showCategory = 5)
p_final <- p + scale_size(range = c(6, 14)) + # increase size
  theme_classic(base_size = 18) +
  theme(plot.title = element_text('test'), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 14), legend.position = "right")
p_final


# plot cluster 4 enrichments

terms_keep <- c('cell surface', 'extracellular matrix', 'endoplasmic reticulum lumen', 'membrane raft', 
                'cell-substrate junction', 'cell-cell junction', 'side of membrane', 'actin cytoskeleton', 'Golgi membrane', 'lytic vacuole')
idx <- GO_results_cluster4@result$Description %in% terms_keep
GO_results_cluster4_selected <- GO_results_cluster4
GO_results_cluster4_selected@result <- GO_results_cluster4@result[idx, ]

p <- dotplot(GO_results_cluster4_selected, showCategory = 5)
p_final <- p + scale_size(range = c(6, 14)) + # increase size
  theme_classic(base_size = 18) +
  theme(plot.title = element_text('test'), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 14), legend.position = "right")
p_final







