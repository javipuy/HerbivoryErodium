##############################################
### E. CICUTARIUM HERBIVORY: FIGURE 4 REDO ###
##############################################

#### Import DESeq2 results ####

setwd("/path/to/project/directory/deseq2/") ## SET CUSTOM PATH TO WORK DIRECTORY

# Unharmed vs mechanical hervibory
data.UNvsMH <- read.table("deseq2_results_ecic-CCvsHC.txt", sep = "\t", header = T, row.names = 1)
# Unharmed vs hervibory with saliva
data.UNvsHS <- read.table("deseq2_results_ecic-CCvsSC.txt", sep = "\t", header = T, row.names = 1)
# Mechanical herbivory vs hervibory with saliva
data.MHvsHS <- read.table("deseq2_results_ecic-HCvsSC.txt", sep = "\t", header = T, row.names = 1)
# Demethylated vs demethylated mechanical hervibory
data.DEvsDM <- read.table("_results_ecic-DEvsDM.txt", sep = "\t", header = T, row.names = 1)
# Demethylated vs demethylated hervibory with saliva
data.DEvsDS <- read.table("_deseq2_results_ecic-DEvsDS.txt", sep = "\t", header = T, row.names = 1)
# Demethylated Mechanical herbivory vs demethylated hervibory with saliva
data.DMvsDS <- read.table("_deseq2_results_ecic-DMvsDS.txt", sep = "\t", header = T, row.names = 1)

#### Generate DESeq2 results for DEvsDM, DEvsDS, and DMvsDS ####

## Uncomment and run this code chunk if DESeq2 results have not bee generated in the session

# # Test Demethylated vs demethylated mechanical hervibory
# res.DEvsDM <- results(dds, contrast = c(0,0,0,1,-1,0))
# data.DEvsDM <- data.frame(res.DEvsDM)
# data.DEvsDM <- data.DEvsDM[with(data.DEvsDM, order(padj)), ]
# write.table(res.DEvsDM, file = "deseq2_results_ecic-DEvsDM.txt", 
#             sep = "\t", row.names = rownames(res.DEvsDM), col.names = colnames(res.DEvsDM), 
#             quote = F, append = F)
# 
# # Test Demethylated vs demethylated hervibory with saliva
# res.DEvsDS <- results(dds, contrast = c(0,0,0,1,0,-1))
# data.DEvsDS <- data.frame(res.DEvsDS)
# data.DEvsDS <- data.DEvsDS[with(data.DEvsDS, order(padj)), ]
# write.table(res.DEvsDS, file = "deseq2_results_ecic-DEvsDS.txt", 
#             sep = "\t", row.names = rownames(res.DEvsDS), col.names = colnames(res.DEvsDS), 
#             quote = F, append = F)
# 
# # Test Demethylated mechanical herbivory vs demethylated hervibory with saliva
# res.DMvsDS <- results(dds, contrast = c(0,0,0,0,1,-1))
# data.DMvsDS <- data.frame(res.DMvsDS)
# data.DMvsDS <- data.DMvsDS[with(data.DMvsDS, order(padj)), ]
# write.table(res.DMvsDS, file = "deseq2_results_ecic-DMvsDS.txt", 
#             sep = "\t", row.names = rownames(res.DMvsDS), col.names = colnames(res.DMvsDS), 
#             quote = F, append = F)

#### DEG lists ####

# Unharmed vs mechanical hervibory
dim(data.UNvsMH[which(data.UNvsMH$padj <= 0.05),]) # 3847 DEGs
dim(data.UNvsMH[which(data.UNvsMH$padj <= 0.05 & data.UNvsMH$log2FoldChange > 0),]) # 2013 Down-regulated by mechanical herbivory
dim(data.UNvsMH[which(data.UNvsMH$padj <= 0.05 & data.UNvsMH$log2FoldChange < 0),]) # 1834 Up-regulated by mechanical herbivory

# Unharmed vs hervibory with saliva
dim(data.UNvsHS[which(data.UNvsHS$padj <= 0.05),]) # 1229 DEGs
dim(data.UNvsHS[which(data.UNvsHS$padj <= 0.05 & data.UNvsHS$log2FoldChange > 0),]) # 853 Down-regulated by hervibory with saliva
dim(data.UNvsHS[which(data.UNvsHS$padj <= 0.05 & data.UNvsHS$log2FoldChange < 0),]) # 376 Up-regulated by hervibory with saliva

# Mechanical herbivory vs hervibory with saliva
dim(data.MHvsHS[which(data.MHvsHS$padj <= 0.05),]) # 3079 DEGs
dim(data.MHvsHS[which(data.MHvsHS$padj <= 0.05 & data.MHvsHS$log2FoldChange > 0),]) # 1554 Up-regulated by mechanical herbivory
dim(data.MHvsHS[which(data.MHvsHS$padj <= 0.05 & data.MHvsHS$log2FoldChange < 0),]) # 1525 Up-regulated by hervibory with saliva

# Demethylated vs demethylated mechanical herbivory
dim(data.DEvsDM[which(data.DEvsDM$padj <= 0.05),]) # 188 DEGs
dim(data.DEvsDM[which(data.DEvsDM$padj <= 0.05 & data.DEvsDM$log2FoldChange > 0),]) # 159 Down-regulated by mechanical herbivory
dim(data.DEvsDM[which(data.DEvsDM$padj <= 0.05 & data.DEvsDM$log2FoldChange < 0),]) # 29 Up-regulated by mechanical herbivory

# Demethylated vs demethylated herbivory with saliva
dim(data.DEvsDS[which(data.DEvsDS$padj <= 0.05),]) # 6 DEGs
dim(data.DEvsDS[which(data.DEvsDS$padj <= 0.05 & data.DEvsDS$log2FoldChange > 0),]) # 3 Down-regulated by herbivory with saliva
dim(data.DEvsDS[which(data.DEvsDS$padj <= 0.05 & data.DEvsDS$log2FoldChange < 0),]) # 3 Up-regulated by herbivory with saliva

# Demethylated mechanical hervibory vs demethylated herbivory with saliva
dim(data.DMvsDS[which(data.DMvsDS$padj <= 0.05),]) # 547 DEGs
dim(data.DMvsDS[which(data.DMvsDS$padj <= 0.05 & data.DMvsDS$log2FoldChange > 0),]) # 371 Up-regulated by mechanical herbivory
dim(data.DMvsDS[which(data.DMvsDS$padj <= 0.05 & data.DMvsDS$log2FoldChange < 0),]) # 176 Up-regulated by herbivory with saliva

# Get DEG lists
UNvsMH.degs <- rownames(data.UNvsMH[which(data.UNvsMH$padj <= 0.05),])
UNvsHS.degs <- rownames(data.UNvsHS[which(data.UNvsHS$padj <= 0.05),])
MHvsHS.degs <- rownames(data.MHvsHS[which(data.MHvsHS$padj <= 0.05),])
UNvsDE.degs <- rownames(data.UNvsDE[which(data.UNvsDE$padj <= 0.05),])
DEvsDM.degs <- rownames(data.DEvsDM[which(data.DEvsDM$padj <= 0.05),])
DEvsDS.degs <- rownames(data.DEvsDS[which(data.DEvsDS$padj <= 0.05),])
DMvsDS.degs <- rownames(data.DMvsDS[which(data.DMvsDS$padj <= 0.05),])

#### Venn Diagrams ####
library(ggvenn) # colours need to be adjusted to alpha = 0.3 in demethylated comparisons
library(gridExtra)

## get pairwise intersections
MHvenn = list(UNvsMH = UNvsMH.degs, DEvsDM = DEvsDM.degs)
MHvenn.plot <- ggvenn(
  MHvenn,
  fill_color = c("#E41A1CFF", "#E41A1C33"),
  stroke_size = 0.5, set_name_size = 4,
  auto_scale = T, show_percentage = F
)

HSvenn = list(UNvsHS = UNvsHS.degs, DEvsDS = DEvsDS.degs)
HSvenn.plot <- ggvenn(
  HSvenn,
  fill_color = c("#377EB8FF", "#377EB833"),
  stroke_size = 0.5, set_name_size = 4,
  auto_scale = T, show_percentage = F
)

MSvenn = list(MHvsHS = MHvsHS.degs, DMvsDS = DMvsDS.degs)
MSvenn.plot <- ggvenn(
  MSvenn,
  fill_color = c("#4DAF4AFF", "#4DAF4A33"),
  stroke_size = 0.5, set_name_size = 4,
  auto_scale = T, show_percentage = F
)

grid.arrange(MHvenn.plot, HSvenn.plot, MSvenn.plot, ncol = 3, nrow = 1)

## get threesome intersections
UNvenn = list(UNvsMH = UNvsMH.degs, UNvsHS = UNvsHS.degs, MHvsHS = MHvsHS.degs)
UNvenn.plot <- ggvenn(
  UNvenn,
  fill_color = c("#E41A1CFF", "#377EB8FF", "#4DAF4AFF"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

DEvenn = list(DEvsDM = DEvsDM.degs, DEvsDS = DEvsDS.degs, DMvsDS = DMvsDS.degs)
DEvenn.plot <- ggvenn(
  DEvenn,
  fill_color = c("#E41A1C33", "#377EB833", "#4DAF4A33"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

grid.arrange(UNvenn.plot, DEvenn.plot, ncol = 2, nrow = 1)


#### Hypergeometric tests ####

#### HYPERGEOMETRIC TEST ####

## Formal hypergeometric test

# autoPhyper function: given 2 DESeq2 output files, give the significance of the overlapping DEGs

autoPhyper <- function(deseqOut1, deseqOut2) {
  deglist1 <- rownames(deseqOut1[which(deseqOut1$padj <= 0.05),])
  group1 <- length(deglist1) # number of DEGs in group 1
  deglist2 <- rownames(deseqOut2[which(deseqOut2$padj <= 0.05),])
  group2 <- length(deglist2) # number of DEGs in group 2
  overlap <- length(intersect(deglist1, deglist2)) # number of overlapping DEGs between groups
  totalgenes <- length(unique(c(rownames(deseqOut1), rownames(deseqOut2)))) # all analyzed genes
  HGtest <- phyper(overlap, group2, totalgenes - group2, group1, lower.tail = FALSE, log.p = FALSE) # Hyper geometric test
  print(paste("Hypergeometric test p-value: ", HGtest)) # print the p-value
}

# UNvsMH, DEvsDM
autoPhyper(data.UNvsMH, data.DEvsDM) # 1.1e-9
# UNvsHS, DEvsDS
autoPhyper(data.UNvsHS, data.DEvsDS) # 0.11
# MHvsHS, DMvsDS
autoPhyper(data.MHvsHS, data.DMvsDS) # 2.7e-81

# UNvsMH, UNvsHS
autoPhyper(data.UNvsMH, data.UNvsHS) # 1.0e-121
# UNvsMH, MHvsHS
autoPhyper(data.UNvsMH, data.MHvsHS) # 0
# UNvsHS, MHvsHS
autoPhyper(data.UNvsHS, data.MHvsHS) # 1.0e-191

# DEvsDM, DEvsDS
autoPhyper(data.DEvsDM, data.DEvsDS) # 1.0e-16
# DEvsDM, DMvsDS
autoPhyper(data.DEvsDM, data.DMvsDS) # 0.001
# DEvsDS, DMvsDS
autoPhyper(data.DEvsDS, data.DMvsDS) # 0.051


#### Volcano Plots ####
library(ggplot2)
library("RColorBrewer")
library(gridExtra)

# UNvsMH
volcDataH <- data.UNvsMH
volcDataH$diffexpr <- "NO"
volcDataH$diffexpr[volcDataH$log2FoldChange > 0 & volcDataH$padj < 0.05] <- "DOWN"
volcDataH$diffexpr[volcDataH$log2FoldChange < 0 & volcDataH$padj < 0.05] <- "UP"
volcPlotH <- ggplot(data = volcDataH, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "UN vs MH") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 FC") + ylab("-log10 FDR") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# UNvsHS
volcDataS <- data.UNvsHS
volcDataS$diffexpr <- "NO"
volcDataS$diffexpr[volcDataS$log2FoldChange > 0 & volcDataS$padj < 0.05] <- "DOWN"
volcDataS$diffexpr[volcDataS$log2FoldChange < 0 & volcDataS$padj < 0.05] <- "UP"
volcPlotS <- ggplot(data = volcDataS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "UN vs HS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 FC") + ylab("-log10 FDR") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# MHvsHS
volcDataHS <- data.MHvsHS
volcDataHS$diffexpr <- "NO"
volcDataHS$diffexpr[volcDataHS$log2FoldChange > 0 & volcDataHS$padj < 0.05] <- "MECHANICAL"
volcDataHS$diffexpr[volcDataHS$log2FoldChange < 0 & volcDataHS$padj < 0.05] <- "SALIVA"
volcPlotHS <- ggplot(data = volcDataHS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "MH vs HS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 FC") + ylab("-log10 FDR") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# all volcano plots in a single figure
grid.arrange(grobs = list(volcPlotH, volcPlotS, volcPlotHS), ncol = 3, nrow = 1) 


#### GO term enrichment analysis plots ####
library(ggplot2)
library(scales)
library(forcats)
library(gridExtra)
# Load these if not loaded already
# library(topGO)
# library(ggplotify)

# load enriched GO terms per comparison
setwd("D:/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_rnaseq_herb/results/topGO/")

goterms <- read.csv("GOtermPlotSummary.csv", sep = ';', header = T)

# generate subsets
goterms.UNvsMHup <- goterms[which(goterms$Comparison == "herbH_up"), ]
goterms.UNvsMHdown <- goterms[which(goterms$Comparison == "herbH_down"), ]
goterms.UNvsHSup <- goterms[which(goterms$Comparison == "herbS_up"), ]
goterms.UNvsHSdown <- goterms[which(goterms$Comparison == "herbS_down"), ]
goterms.MHvsHSup <- goterms[which(goterms$Comparison == "herbHS_H"), ]
goterms.MHvsHSdown <- goterms[which(goterms$Comparison == "herbHS_S"), ]

# generate a funtion to print GO term plots for specific GO term subsets from modified topGO output
goTermPlot <- function(goTermSubset) {
  mainPlot <- 
    ggplot(goTermSubset, 
           aes(x = fct_reorder(Term, weightFDR, .desc = TRUE), 
           y = -log10(weightFDR), 
           size = log2(Significant/Expected), 
           fill = -log10(weightFDR))) +
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'steelblue', high = 'tomato') +
    xlab('') + ylab('Enrichment score') +
    labs(title = paste('Enriched GO terms')) +
    # Theme and general format
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      # Axis format
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      #Legend format
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    coord_flip()
  return(mainPlot)
}
  
goplot.UNvsMHup <- goTermPlot(goterms.UNvsMHup)
goplot.UNvsMHdown <- goTermPlot(goterms.UNvsMHdown)
goplot.UNvsHSup <- goTermPlot(goterms.UNvsHSup)
goplot.UNvsHSdown <- goTermPlot(goterms.UNvsHSdown)
goplot.MHvsHSup <- goTermPlot(goterms.MHvsHSup)
goplot.MHvsHSdown <- goTermPlot(goterms.MHvsHSdown)

library(gridExtra)
goplot.UNvsMHup2 <- ggplotGrob(goplot.UNvsMHup)
goplot.UNvsMHdown2 <- ggplotGrob(goplot.UNvsMHdown)
goplot.UNvsHSup2 <- ggplotGrob(goplot.UNvsHSup)
goplot.UNvsHSdown2 <- ggplotGrob(goplot.UNvsHSdown)
goplot.MHvsHSup2 <- ggplotGrob(goplot.MHvsHSup)
goplot.MHvsHSdown2 <- ggplotGrob(goplot.MHvsHSdown)

# adjust all plot widths to that of goplot.UNvsMHup2
goplot.UNvsMHdown2$widths <- goplot.UNvsMHup2$widths
goplot.UNvsHSup2$widths <- goplot.UNvsMHup2$widths
goplot.UNvsHSdown2$widths <- goplot.UNvsMHup2$widths
goplot.MHvsHSup2$widths <- goplot.UNvsMHup2$widths
goplot.MHvsHSdown2$widths <- goplot.UNvsMHup2$widths

grid.arrange(goplot.UNvsMHup2, goplot.UNvsMHdown2, 
             goplot.UNvsHSup2, goplot.UNvsHSdown2, 
             goplot.MHvsHSup2, goplot.MHvsHSdown2, 
             nrow = 3, ncol = 2) # 20 in x 20 in, rework legends and text manually
