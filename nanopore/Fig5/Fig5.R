library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(ggrastr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd('C:/Users/bacs5867/Documents/PABPN1/scripts/notebooks_publi')

############################################
# Volcano plots

cyto_diff_m6A_rep3_cov <- read.table('Fig5/cyto_diff_m6A_rep3_cov.txt')
cyto_diff_m6A_rep3_cov$log10pval <- -log10(cyto_diff_m6A_rep3_cov$FDR)

volcano <- ggplot(data = cyto_diff_m6A_rep3_cov, aes(delta_prop_m6A, log10pval, col = dir)) + 
  geom_point_rast(size = 2, show.legend = FALSE) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = -0.1, col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = 0.1, col = "gray", linetype = 'dashed') +
  ylab(expression(-log[10]~"(FDR)"))+
  xlab("m6A level delta (shPABPN1 - shCTL)")+
  xlim(-1, 1)+
  ylim(0, 17) +
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "grey"), labels = c("", "", ""), name = '')+
  theme_bw(base_size = 25)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano



cyto_diff_m6A_rep4_cov <- read.table('Fig5/cyto_diff_m6A_rep4_cov.txt')
cyto_diff_m6A_rep4_cov$log10pval <- -log10(cyto_diff_m6A_rep4_cov$FDR)

volcano <- ggplot(data = cyto_diff_m6A_rep4_cov, aes(delta_prop_m6A, log10pval, col = dir)) + 
  geom_point_rast(size = 2, show.legend = FALSE) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = -0.1, col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = 0.1, col = "gray", linetype = 'dashed') +
  ylab(expression(-log[10]~"(FDR)"))+
  xlab("m6A level delta (shPABPN1 - shCTL)")+
  xlim(-1, 1)+
  ylim(0, 17) +
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "grey"), labels = c("", "", ""), name = '')+
  theme_bw(base_size = 25)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano


nucleus_diff_m6A_rep3_cov <- read.table('Fig5/nucleus_diff_m6A_rep3_cov.txt')
nucleus_diff_m6A_rep3_cov$log10pval <- -log10(nucleus_diff_m6A_rep3_cov$FDR)

volcano <- ggplot(data = nucleus_diff_m6A_rep3_cov, aes(delta_prop_m6A, log10pval, col = dir)) + 
  geom_point_rast(size = 2, show.legend = FALSE) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = -0.1, col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = 0.1, col = "gray", linetype = 'dashed') +
  ylab(expression(-log[10]~"(FDR)"))+
  xlab("m6A level delta (shPABPN1 - shCTL)")+
  xlim(-1, 1)+
  ylim(0, 17) +
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "grey"), labels = c("", "", ""), name = '')+
  theme_bw(base_size = 25)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano


###################################################
# GO analysis

cyto_diff_m6A_merge_reps_cov_genes_bg <- read.table('Fig5/cyto_diff_m6A_merge_reps_cov_genes_bg.txt')
cyto_m6A_background <- cyto_diff_m6A_merge_reps_cov_genes_bg$Geneid
cyto_genes_more_m6A <- as.list(read.table('Fig5/cyto_diff_m6A_merge_reps_cov_genes_more_meth.txt')$Geneid)

GO_results_genes_more_m6A_cyto = enrichGO(gene = cyto_genes_more_m6A, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'CC', pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                          universe = cyto_m6A_background)
res_df_genes_more_m6A_cyto <- as.data.frame(GO_results_genes_more_m6A_cyto)

# remove duplicated terms
terms_keep <- c('mitochondrial inner membrane', 'mitochondrial ribosome', 'mitochondrial matrix', 'endoplasmic reticulum lumen', 
                'mitochondrial small ribosomal subunit', 'ribosome', 'pigment granule', 'mitochondrial large ribosomal subunit')
idx <- GO_results_genes_more_m6A_cyto@result$Description %in% terms_keep
GO_results_genes_more_m6A_cyto_selected <- GO_results_genes_more_m6A_cyto
GO_results_genes_more_m6A_cyto_selected@result <- GO_results_genes_more_m6A_cyto@result[idx, ]

p <- dotplot( GO_results_genes_more_m6A_cyto_selected, showCategory = 10)

p_final <- p +
  scale_size(range = c(6, 14)) + # increase size
  theme_classic(base_size = 18) +
  theme(plot.title = element_text('test'), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 14), legend.position = "right")
p_final












