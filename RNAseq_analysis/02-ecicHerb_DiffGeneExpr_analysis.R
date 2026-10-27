#### E. CICUTARIUM HERBIVORY ####

setwd("/path/to/project/directory/deseq2") ## SET CUSTOM PATH TO WORK DIRECTORY

library(DESeq2)
counts <- read.csv("../featureCounts/read-counts-greater100.csv", sep = ";", header = T, row.names = 1)

library(stringr)
samples <- colnames(counts)
herbivory <- as.factor(str_sub(samples, -3, -3))
azacytidine <- as.factor(str_sub(samples, -1, -1))
individual <- as.factor(str_sub(samples, 5, 6))
group <- as.factor(str_sub(samples, -3, -1))
herbivory <- relevel(herbivory, "C")
azacytidine <- relevel(azacytidine, "C")
group <- relevel(group, "C_C")

table <- data.frame(sampleName = samples, fileName = samples, condition = herbivory)
table$azacytidine <- azacytidine

#### For herbivory and azacitidine general effects ####

expDesign <- model.matrix(~ herbivory + azacytidine + herbivory:azacytidine) 
ddsMatrix <- DESeqDataSetFromMatrix(countData = counts, colData = table, design = expDesign)

dds <- ddsMatrix[rowSums(counts(ddsMatrix)) > 1, ]
rld <- rlog(dds, blind = FALSE)
dds <- DESeq(ddsMatrix)
normCounts <- counts(dds, normalized = TRUE)
write.table(normCounts, file = "ecic-herb_normcounts_deseq2.txt", 
            sep = "\t", row.names = rownames(normCounts), col.names = colnames(normCounts), 
            quote = F, append = F)

par(mfcol = c(2,1))
boxplot(log(counts,2), col = as.factor(group)) # Raw counts
boxplot(log(normCounts, 2), col = as.factor(group)) # Norm counts

## Explanation of the model factor and levels
# Samples can be a combination of the factor *herbivory*, with levels: 
# -C(ontrol, undamaged)
# -H (Mechanical damage)
# -S (Herbivory simulated with Saliva)
# and the factor *azacitidine*, with levels:
# -C(ontrol, not demethylated)
# -A(zacytidine, demethylated)
# Thus, the combined factor levels present two columns: 
# one specifying either C, M, or H, from herbivory factor;
# another specifying C, or A, from azacitidine factor.
# When fitting the model, the reference treatment would be CC (Undamaged Controls)
# When specifying which comparison frm the DGEA are we obtaining results from,
# we use a vector specifying which levels are we comparing.
# The vector (specified in the contrast option), will show six positions:
# -The first one is the intercept.
# -second to sixth positions are levels HC, SC, CA, HA, SA.
# The positions may show three values:
# -zero (0): this level is not used in the comparison
# -plus one (1): this is the treatment, placed in the numerator when calculating fold changes or stats.
# -minus one (-1): this is the control or reference set of samples, placed in the denominator
# This way, if we want to compare unharmed controls (C_C) with unharmed demethylated (C_A) samples,
# we would specify in the contrast vector the following:
# results(dds, contrast = c(0,0,0,-1,0,0)) 
# The CC treatment would be the reference, thus we do not modify the vector, but we put -1 in the position which C_A treatment is (the fourth)
# If we wanted to compare the mechanichally damaged controls (HC) with the saliva induced demethylated samples (SA),
# we would specify them as follows: 
# results(dds, contrast = c(0,-1,0,0,0,1))
# This time, we give both the treatment (1) and the reference (-1) set of samples.
# To test all control vs all demethylated samples, this would be the command:  
# results(dds, contrast = c(0,-1,-1,1,1,1))
# With all the controls with -1, and all the demethylated sample with 1. We leave the first position (intercept) as zero.

# Unharmed vs mechanical hervibory
res.CCvsHC <- results(dds, contrast = c(0,-1,0,0,0,0))
data.CCvsHC <- data.frame(res.CCvsHC)
data.CCvsHC <- data.CCvsHC[with(data.CCvsHC, order(padj)), ]
write.table(res.CCvsHC, file = "deseq2_results_ecic-CCvsHC.txt", 
            sep = "\t", row.names = rownames(res.CCvsHC), col.names = colnames(res.CCvsHC), 
            quote = F, append = F)
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),]) # 3847 DEGs
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange > 0),]) # 2013 Down-regulated herbivoryH
dim(data.CCvsHC[which(data.CCvsHC$padj <= 0.05 & data.CCvsHC$log2FoldChange < 0),]) # 1834 Up-regulated herbivoryH

# Unharmed vs hervibory with saliva
res.CCvsSC <- results(dds, contrast = c(0,0,-1,0,0,0))
data.CCvsSC <- data.frame(res.CCvsSC)
data.CCvsSC <- data.CCvsSC[with(data.CCvsSC, order(padj)), ]
write.table(res.CCvsSC, file = "deseq2_results_ecic-CCvsSC.txt", 
            sep = "\t", row.names = rownames(res.CCvsSC), col.names = colnames(res.CCvsSC), 
            quote = F, append = F)
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]) # 1229 DEGs
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange > 0),]) # 853 Down-regulated herbivoryS
dim(data.CCvsSC[which(data.CCvsSC$padj <= 0.05 & data.CCvsSC$log2FoldChange < 0),]) # 376 Up-regulated herbivoryS

# Mechanical hervibory vs hervibory with saliva
res.HCvsSC <- results(dds, contrast = c(0,1,-1,0,0,0))
data.HCvsSC <- data.frame(res.HCvsSC)
data.HCvsSC <- data.HCvsSC[with(data.HCvsSC, order(padj)), ]
write.table(res.HCvsSC, file = "deseq2_results_ecic-HCvsSC.txt", 
            sep = "\t", row.names = rownames(res.HCvsSC), col.names = colnames(res.HCvsSC), 
            quote = F, append = F)
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]) # 3079 DEGs
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange > 0),]) # 1554 Up-regulated herbivoryH
dim(data.HCvsSC[which(data.HCvsSC$padj <= 0.05 & data.HCvsSC$log2FoldChange < 0),]) # 1525 Up-regulated herbivoryS


# Unharmed Control vs Unharmed Demethylated
res.CCvsCA <- results(dds, contrast = c(0,0,0,-1,0,0))
data.CCvsCA <- data.frame(res.CCvsCA)
data.CCvsCA <- data.CCvsCA[with(data.CCvsCA, order(padj)), ]
write.table(res.CCvsCA, file = "deseq2_results_ecic-CCvsCA.txt", 
            sep = "\t", row.names = rownames(res.CCvsCA), col.names = colnames(res.CCvsCA), 
            quote = F, append = F)
dim(data.CCvsCA[which(data.CCvsCA$padj <= 0.05),]) # 40 DEGs
dim(data.CCvsCA[which(data.CCvsCA$padj <= 0.05 & data.CCvsCA$log2FoldChange > 0),]) # 14 Down-regulated Azacytidine
dim(data.CCvsCA[which(data.CCvsCA$padj <= 0.05 & data.CCvsCA$log2FoldChange < 0),]) # 26 Up-regulated Azacytidine


## Check expression of top DEG for a contrast

# black = CC, red = CA, purple = SC, cyan = SA, blue = HC, green = HA
# Positive logFC
DEGs.p <- data.CvsA[which(data.CvsA$padj <= 0.05 & data.CvsA$log2FoldChange > 0),]
topDEG.p <- rownames(DEGs.p[1,])
topDEG.p.normCounts <- normCounts[which(rownames(normCounts) == topDEG.p),]
barplot(topDEG.p.normCounts, col = as.factor(group))

# Negative logFC
DEGs.m <- data.CvsA[which(data.CvsA$padj <= 0.05 & data.CvsA$log2FoldChange < 0),]
topDEG.m <- rownames(DEGs.m[1,])
topDEG.m.normCounts <- normCounts[which(rownames(normCounts) == topDEG.m),]
barplot(topDEG.m.normCounts, col = as.factor(group))


#### PLOTS ####

### NMDS ###

sampleColor <- c("white","red","white","red","white","red",
                 "red","white","red","white","orange","orange",
                 "orange","brown","brown","orange","brown","orange",
                 "brown","orange","brown","brown","red","green",
                 "green","green","green","green","green","darkgreen",
                 "darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","white") 

library(edgeR)
plotMDS(dds, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, #labels = colnames(dds), 
        pch = 21,
        bg = sampleColor)
dev.copy(pdf, 'ecic-herb_deseq2_pca.pdf') # C_C white,C_A red,H_C green,H_A darkgreen,S_C orange,S_A brown
dev.off()

# load data to the R session if not already loaded
# normCounts <- read.table("ecic-herb_normcounts_deseq2.txt", sep = "\t", header = T, row.names = 1)
# data.CCvsHC <- read.table("deseq2_results_ecic-CCvsHC.txt", sep = "\t", header = T, row.names = 1)
# data.CCvsSC <- read.table("deseq2_results_ecic-CCvsSC.txt", sep = "\t", header = T, row.names = 1)
# data.HCvsSC <- read.table("deseq2_results_ecic-HCvsSC.txt", sep = "\t", header = T, row.names = 1)
# data.CCvsCA <- read.table("deseq2_results_ecic-CCvsCA.txt", sep = "\t", header = T, row.names = 1)

# PCA for CCvsHC DEGs
DEGs.CCvsHC <- rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])
DEGs.normCounts.CCvsHC <- normCounts[which(DEGs.CCvsHC %in% rownames(normCounts)), ]
# PCA for CCvsSC DEGs
DEGs.CCvsSC <- rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),])
DEGs.normCounts.CCvsSC <- normCounts[which(DEGs.CCvsSC %in% rownames(normCounts)), ]
# PCA for HCvsSC DEGs
DEGs.HCvsSC <- rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),])
DEGs.normCounts.HCvsSC <- normCounts[which(DEGs.HCvsSC %in% rownames(normCounts)), ]
# PCA for CCvsCA DEGs
DEGs.CCvsCA <- rownames(data.CCvsCA[which(data.CCvsCA$padj <= 0.05),])
DEGs.normCounts.CCvsCA <- normCounts[which(DEGs.CCvsCA %in% rownames(normCounts)), ]

par(mfcol = c(2,2))
plotMDS(DEGs.normCounts.CCvsHC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)
plotMDS(DEGs.normCounts.CCvsSC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)
plotMDS(DEGs.normCounts.HCvsSC, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)
plotMDS(DEGs.normCounts.CCvsCA, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)

# PCA for all DEGs
DEGs.all <- unique(c(DEGs.CCvsHC, DEGs.CCvsSC, DEGs.HCvsSC, DEGs.CCvsCA))
DEGs.normCounts.all <- normCounts[which(DEGs.all %in% rownames(normCounts)), ]
plotMDS(DEGs.normCounts.all, gene.selection = "common", top = 6419, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)


### Number of Up-regulated and Down-regulated DEGs with ggplot2 ###

# HerbivoryH
geneExpCH <- cbind(c(1834, 2013))
colnames(geneExpCH) <- "herbivoryH"
rownames(geneExpCH) <- c("Up-regulated", "Down-regulated")

library(reshape2)
geneExpCH.melt <- melt(geneExpCH)

library(ggplot2)
plotCH.2 <- ggplot(data = geneExpCH.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Herbivory") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# HerbivoryS
geneExpCS <- cbind(c(376, 853))
colnames(geneExpCS) <- "herbivoryS"
rownames(geneExpCS) <- c("Up-regulated", "Down-regulated")

library(reshape2)
geneExpCS.melt <- melt(geneExpCS)

library(ggplot2)
plotCS.2 <- ggplot(data = geneExpCS.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Saliva") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Azacytidine
geneExpAC <- cbind(c(26, 14))
colnames(geneExpAC) <- "azacytidine"
rownames(geneExpAC) <- c("Up-regulated", "Down-regulated")

library(reshape2)
geneExpAC.melt <- melt(geneExpAC)

library(ggplot2)
plotAC.2 <- ggplot(data = geneExpAC.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Azacitidine") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Within Herbivory
geneExpHS <- cbind(c(1554, 1525))
colnames(geneExpHS) <- "withinHerbivory"
rownames(geneExpHS) <- c("HerbivoryH.Up", "HerbivoryS.Up")

library(reshape2)
geneExpHS.melt <- melt(geneExpHS)

library(ggplot2)
plotHS.2 <- ggplot(data = geneExpHS.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Herbivory H vs S") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")


library(gridExtra)
grid.arrange(plotCH.2, plotCS.2, plotHS.2, plotAC.2, ncol = 2, nrow = 2) # plotAA.2 left out


### HEATMAP ###
# DEGs overlapping the three contrasts
# Some columns added in Excel

# Set DEG lists
DEGs.CCvsHC <- data.CCvsHC[which(data.CCvsHC$padj <= 0.05),]
DEGs.CCvsSC <- data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]
DEGs.HCvsSC <- data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]
DEGs.CCvsCA <- data.CCvsCA[which(data.CCvsCA$padj <= 0.05),]

# Set DEGs normalized counts
normCounts.CCvsHC <- normCounts[which(rownames(normCounts) %in% rownames(DEGs.CCvsHC)), ]
normCounts.CCvsSC <- normCounts[which(rownames(normCounts) %in% rownames(DEGs.CCvsSC)), ]
normCounts.HCvsSC <- normCounts[which(rownames(normCounts) %in% rownames(DEGs.HCvsSC)), ]
normCounts.CCvsCA <- normCounts[which(rownames(normCounts) %in% rownames(DEGs.CCvsCA)), ]

# Merge DEG normCounts and remove duplicated rows
degs <- rbind(normCounts.CCvsHC, normCounts.CCvsSC, normCounts.HCvsSC, normCounts.CCvsCA)
degs <- degs[!duplicated(degs), ]

library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)

samples.h <- colnames(degs)
genes.h <- rownames(degs)
group.h <- as.factor(paste0(herbivory, azacytidine))
levels(group.h) <- c("black","red","green","darkgreen","orange","brown")
# C_C black,C_A red,H_C green,H_A darkgreen,S_C orange,S_A brown

library(gplots)
heatmap.2(as.matrix(scale(degs)), trace = "none", col = coul, scale = "row", labRow = FALSE, colCol = as.vector(group.h))


### Volcano Plots ###

library(ggplot2)
library("RColorBrewer")

# HerbivoryH
volcDataH <- read.table('deseq2_results_ecic-CCvsHC.txt', sep = '\t', header = TRUE)
volcDataH$diffexpr <- "NO"
volcDataH$diffexpr[volcDataH$log2FoldChange > 0 & volcDataH$padj < 0.05] <- "DOWN"
volcDataH$diffexpr[volcDataH$log2FoldChange < 0 & volcDataH$padj < 0.05] <- "UP"
volcPlotH <- ggplot(data = volcDataH, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryH") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# HerbivoryS
volcDataS <- read.table('deseq2_results_ecic-CCvsSC.txt', sep = '\t', header = TRUE)
volcDataS$diffexpr <- "NO"
volcDataS$diffexpr[volcDataS$log2FoldChange > 0 & volcDataS$padj < 0.05] <- "DOWN"
volcDataS$diffexpr[volcDataS$log2FoldChange < 0 & volcDataS$padj < 0.05] <- "UP"
volcPlotS <- ggplot(data = volcDataS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))

# Azacytidine
volcDataA <- read.table('deseq2_results_ecic-CCvsCA.txt', sep = '\t', header = TRUE)
volcDataA$diffexpr <- "NO"
volcDataA$diffexpr[volcDataA$log2FoldChange > 0 & volcDataA$padj < 0.05] <- "DOWN"
volcDataA$diffexpr[volcDataA$log2FoldChange < 0 & volcDataA$padj < 0.05] <- "UP"
volcPlotA <- ggplot(data = volcDataA, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs Azacytidine") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-15, 15)) + scale_y_continuous(limits = c(0, 7.5))

# Within Herbivory
volcDataHS <- read.table('deseq2_results_ecic-HCvsSC.txt', sep = '\t', header = TRUE)
volcDataHS$diffexpr <- "NO"
volcDataHS$diffexpr[volcDataHS$log2FoldChange > 0 & volcDataHS$padj < 0.05] <- "HERBIVORY"
volcDataHS$diffexpr[volcDataHS$log2FoldChange < 0 & volcDataHS$padj < 0.05] <- "SALIVA"
volcPlotHS <- ggplot(data = volcDataHS, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "HerbivoryH vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 25))


library(gridExtra)
grid.arrange(volcPlotH, volcPlotS, volcPlotHS, volcPlotA, ncol = 2, nrow = 2) # volcPlotAA left out

### MA PLOTS ###

MAPlotH <- ggplot(data = volcDataH, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryH") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 fold change") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

MAPlotS <- ggplot(data = volcDataS, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 fold change") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

MAPlotA <- ggplot(data = volcDataA, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs Azacytidine") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 fold change") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

MAPlotHS <- ggplot(data = volcDataHS, aes(x = log2(baseMean), y = log2FoldChange, colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "HerbivoryH vs HerbivoryS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 CPM") + ylab("log2 fold change") +
  scale_colour_manual(values = c("steelblue", "grey", "red"))

library(gridExtra)
grid.arrange(MAPlotH, MAPlotS, MAPlotHS, MAPlotA, ncol = 2, nrow = 2)

#### HYPERGEOMETRIC TEST ####

## Intersections
# DEGs between CCvsCA and CCvsHC: 6
length(intersect(as.vector(rownames(data.CCvsCA[which(data.CCvsCA$padj <= 0.05),])), as.vector(rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),]))))
# DEGs between CCvsCA and CCvsSC: 8
length(intersect(as.vector(rownames(data.CCvsCA[which(data.CCvsCA$padj <= 0.05),])), as.vector(rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]))))
# DEGs between CCvsCA and HCvsSC: 1
length(intersect(as.vector(rownames(data.CCvsCA[which(data.CCvsCA$padj <= 0.05),])), as.vector(rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]))))
# DEGs between CCvsHC and CCvsSC: 331
length(intersect(as.vector(rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])), as.vector(rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),]))))
# DEGs between CCvsHC and HCvsSC: 1084
length(intersect(as.vector(rownames(data.CCvsHC[which(data.CCvsHC$padj <= 0.05),])), as.vector(rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]))))
# DEGs between CCvsSC and HCvsSC: 375
length(intersect(as.vector(rownames(data.CCvsSC[which(data.CCvsSC$padj <= 0.05),])), as.vector(rownames(data.HCvsSC[which(data.HCvsSC$padj <= 0.05),]))))

## Formal hypergeometric test
# CCvsCA U CCvsHC
phyper(5, 40, 62477, 3847, lower.tail = FALSE, log.p = FALSE) # P-value < 0.05 (0.03433348) *
# CCvsCA U CCvsSC
phyper(7, 40, 62477, 1229, lower.tail = FALSE, log.p = FALSE) # P-value < 0.001 (9.595232e-07) ***
# CCvsCA U HCvsSC
phyper(0, 40, 62477, 3079, lower.tail = FALSE, log.p = FALSE) # P-value > 0.05 (0.8674556)
# CCvsHC U CCvsSC
phyper(330, 3847, 58670, 1229, lower.tail = FALSE, log.p = FALSE) # P-value < 0.001 (6.34611e-121) ***
# CCvsHC U HCvsSC
phyper(1083, 3847, 58670, 3079, lower.tail = FALSE, log.p = FALSE) # P-value < 0.001 (0) ***
# CCvsSC U HCvsSC
phyper(374, 1229, 61288, 3079, lower.tail = FALSE, log.p = FALSE) # P-value < 0.001 (9.556417e-191) ***

#q = size of overlap - 1
#m = number of DEGs in experiment comparison 1
#n = total number of genes (62517) - m
#k = number of DEGs in comparison experiment 2


#### PLOTS FOR GOIs ####

# Annotated DEGs in the Control vs Azacitidine comparison
annotationsDEG <- c("g25155.t1", "g22465.t1", "g81973.t1", "g64755.t1",	
                    "g169477.t1",	"g121748.t1",	"g15886.t1", "g141517.t1", 
                    "g148881.t1", "g61285.t1",	"g178124.t1",	"g60462.t1", 
                    "g88692.t1",	"g95990.t1")
degs.annot <- degs[, annotationsDEG]
degs.annot$Herbivory <- degs$Herbivory
degs.annot$Azacitidine <- degs$Azacitidine

library(reshape2)
degs.annot.melt <- melt(degs.annot)
library(ggplot2)
degs.annot.plot <- ggplot(degs.annot.melt, aes(x = Herbivory, y = log(value), color = Azacitidine)) + 
  #geom_line(aes(group = treatment)) + 
  geom_boxplot() + #scale_x_discrete(expand = c(0.2, 0.2)) + 
  facet_wrap(~ variable, scale = "free")

