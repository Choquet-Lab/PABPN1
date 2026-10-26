library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(ggrastr)

setwd('C:/Users/bacs5867/Documents/PABPN1/scripts/notebooks_publi')

ensembl_id_gene_name <- read.table('C:/Users/bacs5867/Documents/PABPN1/analyse/hg38_genome_ensembl_114/Homo_sapiens.GRCh38.114_ensembl_gene_names.txt')

protein_coding_genes <- read.table('Fig4/biomart_ensembl_mrna.txt', header = T)
protein_coding_genes$mRNA = T

lncrna_genes <- read.table('Fig4/biomart_ensembl_lncrna.txt', header = T)
lncrna_genes$lncrna = T

tpm <- read.table('Fig4/tpm.txt')[,c(1,10,11)]


#########################################
# cytoplasm diff expression

counts_cyto <- read.table('Fig4/cyto_readcounts.txt', header = T)
counts_cyto <- merge(counts_cyto, tpm, on = 'Geneid')
counts_cyto <- counts_cyto[counts_cyto$expression_cyto == 'high',] # filter out genes with less than 1 tpm in all cyto samples of one condition
row.names(counts_cyto) <- counts_cyto$Geneid
counts_cyto <- counts_cyto[,c(2,3,4,5)]
coldata_cyto <- read.table('Fig4/cyto_coldata.txt', header = T)

dds_cyto <- DESeqDataSetFromMatrix(countData = counts_cyto, colData = coldata_cyto, design = ~ condition)
deseq_cyto <- DESeq(dds_cyto)
res_cyto <- results(deseq_cyto, alpha = 0.05)
res_cyto
summary(res_cyto)

resDF_cyto <- as.data.frame(res_cyto)
resDF_cyto <- resDF_cyto[complete.cases(resDF_cyto),] # remove rows with Na
resDF_cyto$Geneid = rownames(resDF_cyto)
resDF_cyto <- merge(resDF_cyto, ensembl_id_gene_name, by = 'Geneid', all.x = T)
resDF_cyto <- resDF_cyto[order(resDF_cyto$padj),]

resDF_cyto$gene_name[is.na(resDF_cyto$gene_name)] <- resDF_cyto$Geneid[is.na(resDF_cyto$gene_name)] # get ensembl id where no gene name from gtf
resDF_cyto$log10pval <- -log10(resDF_cyto$padj)

resDF_cyto$diff_expression <- 'NO'
resDF_cyto$diff_expression[resDF_cyto$log2FoldChange>1 & resDF_cyto$padj<0.05] <- 'UP'
resDF_cyto$diff_expression[resDF_cyto$log2FoldChange<(-1) & resDF_cyto$padj<0.05] <- 'DOWN'

resDF_cyto <- merge(resDF_cyto, protein_coding_genes, by = 'Geneid', all.x = T)
resDF_cyto <- merge(resDF_cyto, lncrna_genes, by = 'Geneid', all.x = T)
resDF_cyto$type <- 'other'
resDF_cyto$type[resDF_cyto$lncrna == T] <- 'lncrna'
resDF_cyto$type[resDF_cyto$mRNA == T] <- 'mrna'
resDF_cyto <- resDF_cyto[,-c(11,12)]

resDF_cyto <- resDF_cyto[order(resDF_cyto$padj),]
write.table(resDF_cyto, 'Fig4/resDF_cyto.txt', sep = '\t')

##############################################
# nucleus diff expression

counts_nucleus <- read.table('Fig4/nucleus_readcounts.txt', header = T)
counts_nucleus <- merge(counts_nucleus, tpm, on = 'Geneid')
counts_nucleus <- counts_nucleus[counts_nucleus$expression_nucleus == 'high',] # filter out genes with less than 1 tpm in all nucleus samples of one condition
row.names(counts_nucleus) <- counts_nucleus$Geneid
counts_nucleus <- counts_nucleus[,c(2,3,4,5)]
coldata_nucleus <- read.table('Fig4/nucleus_coldata.txt', header = T)

dds_nucleus <- DESeqDataSetFromMatrix(countData = counts_nucleus, colData = coldata_nucleus, design = ~ condition)
deseq_nucleus <- DESeq(dds_nucleus)
res_nucleus <- results(deseq_nucleus, alpha = 0.05)
res_nucleus
summary(res_nucleus)

resDF_nucleus <- as.data.frame(res_nucleus)
resDF_nucleus <- resDF_nucleus[complete.cases(resDF_nucleus),] # remove rows with Na
resDF_nucleus$Geneid = rownames(resDF_nucleus)
resDF_nucleus <- merge(resDF_nucleus, ensembl_id_gene_name, by = 'Geneid', all.x = T)
resDF_nucleus <- resDF_nucleus[order(resDF_nucleus$padj),]

resDF_nucleus$gene_name[is.na(resDF_nucleus$gene_name)] <- resDF_nucleus$Geneid[is.na(resDF_nucleus$gene_name)] # get ensembl id where no gene name from gtf
resDF_nucleus$log10pval <- -log10(resDF_nucleus$padj)

resDF_nucleus$diff_expression <- 'NO'
resDF_nucleus$diff_expression[resDF_nucleus$log2FoldChange>1 & resDF_nucleus$padj<0.05] <- 'UP'
resDF_nucleus$diff_expression[resDF_nucleus$log2FoldChange<(-1) & resDF_nucleus$padj<0.05] <- 'DOWN'

resDF_nucleus <- merge(resDF_nucleus, protein_coding_genes, by = 'Geneid', all.x = T)
resDF_nucleus <- merge(resDF_nucleus, lncrna_genes, by = 'Geneid', all.x = T)
resDF_nucleus$type <- 'other'
resDF_nucleus$type[resDF_nucleus$lncrna == T] <- 'lncrna'
resDF_nucleus$type[resDF_nucleus$mRNA == T] <- 'mrna'
resDF_nucleus <- resDF_nucleus[,-c(11,12)]

resDF_nucleus <- resDF_nucleus[order(resDF_nucleus$padj),]
write.table(resDF_nucleus, 'Fig4/resDF_nucleus.txt', sep = '\t')


##############################################
# ratio cyto/nucleus

counts <- read.table('Fig4/readcounts.txt', header = T)
row.names(counts) <- counts$Geneid
counts <- counts[,c(2,3,4,5,6,7,8,9)]
coldata <- read.table('Fig4/coldata.txt', header = T)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition + compartment + condition:compartment)
deseq <- DESeq(dds)
resultsNames(deseq)
res <- results(deseq, alpha = 0.05, name="conditionshPABPN1.compartmentnucleus")
res
summary(res)
resDF <- as.data.frame(res)

resDF$Geneid = rownames(resDF)
resDF <- merge(resDF, ensembl_id_gene_name, by = 'Geneid', all.x = T)
resDF$gene_name[is.na(resDF$gene_name)] <- resDF$Geneid[is.na(resDF$gene_name)] # get ensembl id where no gene name from gtf
resDF <- merge(resDF, protein_coding_genes, by = 'Geneid', all.x = T)
resDF <- merge(resDF, lncrna_genes, by = 'Geneid', all.x = T)
resDF <- resDF[order(resDF$padj),]

resDF$type <- 'other'
resDF$type[resDF$lncrna == T] <- 'lncRNA'
resDF$type[resDF$mRNA == T] <- 'mRNA'

write.table(resDF, 'Fig4/resDF_all.txt', sep = '\t')
