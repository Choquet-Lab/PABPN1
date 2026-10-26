library(APALORD)
library(GenomeInfoDb)
library(txdbmaker)
library(GenomicAlignments)
library(data.table)
library(ggplot2)
library(ggrastr)
library(RColorBrewer)

setwd('C:/Users/bacs5867/Documents/PABPN1/scripts/notebooks_publi')

##########################################


path <- 'C:/Users/bacs5867/Documents/PABPN1/isoquant/dorado_new'

fasta_ref <- 'C:/Users/bacs5867/Documents/PABPN1/analyse/Homo_sapiens.GRCh38.dna.primary_assembly.and.sacCer3_ENO2.fa'
gtf_ref <- 'C:/Users/bacs5867/Documents/PABPN1/analyse/Homo_sapiens.GRCh38.110.gtf'
gene_reference <- load_gtf(gtf_ref)


# cytoplasm
shCTL_cyto <- c(paste0(path, "/shCTL_cyto_rep2"), paste0(path, "/shCTL_cyto_rep3"), paste0(path, "/shNT_cyto"))
shPABPN1_cyto <- c(paste0(path, "/shPAB_cyto_rep2"), paste0(path, "/shPAB_cyto_rep3"), paste0(path, "/shPAB_cyto"))

reads <- load_samples(shCTL_cyto, shPABPN1_cyto, group1 = "shCTL_cyto", group2 = 'shPABPN1_cyto')

# gene level analysis
APA_data <- APA_profile(gene_reference, reads, control = "shCTL_cyto", experimental = "shPABPN1_cyto", direct_RNA = TRUE)
APA_gene_table <- APA_plot(APA_data)
APA_gene_table$log10padj = -log10(APA_gene_table$P_adj)
APA_gene_table <- APA_gene_table[order(APA_gene_table$P_adj),]
write.table(APA_gene_table, file = 'Fig2/APA_gene_table_cyto.txt', sep = "\t")

# site level analysis
PAU_data <- PAU_by_sample(gene_reference, reads, direct_RNA = TRUE)
PAU_test_data <- PAU_test(PAU_data, reads, P_cutoff = 0.05)
write.table(PAU_test_data, file = 'Fig2/cyto_PAU_stats.txt', sep = "\t")



# nucleus
shCTL_nucleus <- c(paste0(path, "/shCTL_noyau_rep2"), paste0(path, "/shCTL_noyau_rep3"), paste0(path, "/shNT_noyau"))
shPABPN1_nucleus <- c(paste0(path, "/shPAB_noyau_rep2"), paste0(path, "/shPAB_noyau_rep3"), paste0(path, "/shPAB_noyau"))

reads <- load_samples(shCTL_nucleus, shPABPN1_nucleus, group1 = "shCTL_nucleus", group2 = 'shPABPN1_nucleus')

# gene level analysis
APA_data <- APA_profile(gene_reference, reads, control = "shCTL_nucleus", experimental = "shPABPN1_nucleus", direct_RNA = TRUE)
APA_gene_table <- APA_plot(APA_data)
APA_gene_table$log10padj = -log10(APA_gene_table$P_adj)
APA_gene_table <- APA_gene_table[order(APA_gene_table$P_adj),]
write.table(APA_gene_table, file = 'Fig2/APA_gene_table_nucleus.txt', sep = "\t")

# site level analysis
PAU_data <- PAU_by_sample(gene_reference, reads, direct_RNA = TRUE)
PAU_test_data <- PAU_test(PAU_data, reads, P_cutoff = 0.05)
write.table(PAU_test_data, file = 'Fig2/nucleus_PAU_stats.txt', sep = "\t")



















