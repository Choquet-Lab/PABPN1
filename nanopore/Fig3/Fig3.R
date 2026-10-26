library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(ggrastr)

setwd('C:/Users/bacs5867/Documents/PABPN1/scripts/notebooks_publi')


##### Long reads

nucleus_introns_volcano <- read.table('Fig3/nucleus_introns_volcano.txt', sep = '\t', header = T)

volcano <- ggplot(data = nucleus_introns_volcano, aes(delta, log10pvalue, col = dir)) + 
  geom_point_rast(size = 1, raster.dpi = 300) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  ggtitle('Nucleus')+
  ylab(expression(-log[10]~"(p-value)"))+
  xlab('IR ratio delta')+
  xlim(-1, 1)+
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "grey"))+
  theme_bw(base_size = 20)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano


cyto_introns_volcano <- read.table('Fig3/cyto_introns_volcano.txt', sep = '\t', header = T)

volcano <- ggplot(data = cyto_introns_volcano, aes(delta, log10pvalue, col = dir)) + 
  geom_point_rast(size = 1, raster.dpi = 300) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  ggtitle('Cytoplasm')+
  ylab(expression(-log[10]~"(p-value)"))+
  xlab('IR ratio delta')+
  xlim(-1, 1)+
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "grey"))+
  theme_bw(base_size = 20)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano



##### Short reads


shPABPN1_vs_shNT_nucleus_volcano <- read.table('Fig3/shPABPN1_vs_shNT_nucleus_volcano.txt', sep = '\t', header = T)

volcano <- ggplot(data = shPABPN1_vs_shNT_nucleus_volcano, aes(delta_IR, log10pvalue, col = dir)) + 
  geom_point_rast(size = 1, raster.dpi = 300) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  ggtitle('Nucleus')+
  ylab(expression(-log[10]~"(p-value)"))+
  xlab('IR ratio delta')+
  xlim(-1, 1)+
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "grey"))+
  theme_bw(base_size = 20)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano


shPABPN1_vs_shNT_cyto_volcano <- read.table('Fig3/shPABPN1_vs_shNT_cyto_volcano.txt', sep = '\t', header = T)

volcano <- ggplot(data = shPABPN1_vs_shNT_cyto_volcano, aes(delta_IR, log10pvalue, col = dir)) + 
  geom_point_rast(size = 1, raster.dpi = 300) + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  ggtitle('Cytoplasm')+
  ylab(expression(-log[10]~"(p-value)"))+
  xlab('IR ratio delta')+
  xlim(-1, 1)+
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "grey"))+
  theme_bw(base_size = 20)+
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3))
volcano
