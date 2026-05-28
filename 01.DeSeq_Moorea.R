######################################Set up######################################
### Load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("genefilter")
library(genefilter)

BiocManager::install("DEGreport")
library(DEGreport)

BiocManager::install("vsn")
library(vsn)

BiocManager::install("apeglm")
library(apeglm)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

install.packages("devtools")
library(devtools)

install_github("stephens999/ashr")
library(ashr)

install.packages("factoextra")
library(factoextra)

install.packages("tidyverse")
library(tidyverse)

library(stringr)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(dendextend)
library(ggplot2)
library(gplots)
library(gridExtra)
library(ggrepel)

### Clean up the workspace
rm(list = ls()) 

### Set directory
wd <- "/Users/marleenklann/Desktop/manini/developmental_transcriptomics/01.DeSeq/zoe/" 
setwd(wd)

### Read file
data <- read.table("../../00.degnorm_counts/zoe/adjusted_read_counts_all_samples.txt", header = TRUE, sep = "\t", dec = ".") 

### Format data frame
rownames(data) <- data$target_id
data <- data[, -1]

### Remove weird samples
data <- data %>% dplyr::select(-BR_JUV_ZC0113, -MAN_D8_ZC0189, -MAN_JUV_ZC0197, -SB_JUV_ZC0168)  

### Get counts
countData <- data %>% round()
write.table(countData, "./results/raw.counts.txt", sep = "\t")

### Read metadata
sampleInfo <- read.table("./metadata.txt", header = TRUE, sep = "\t")
rownames(sampleInfo) <- sampleInfo$sample
sampleInfo$sample <- NULL

rows_to_remove <- c("BR_JUV_ZC0113", "MAN_D8_ZC0189", "MAN_JUV_ZC0197", "SB_JUV_ZC0168")

# Remove rows based on the condition
sampleInfo <- subset(sampleInfo, !(row.names(sampleInfo) %in% rows_to_remove))

### Verify metadata
all(rownames(sampleInfo) == colnames(countData))

### Add new column
sampleInfo$condition <- paste(sampleInfo$site, sampleInfo$day, sep = "_")

## In case you need to select subset of samples
#juv_cols <- grep("JUV", names(countData), value = TRUE)
#juv <- countData[juv_cols]
#juv_sampleInfo <- sampleInfo[sampleInfo$day == "Juvenile", ]


######################################1. QC Analyses######################################
### Check how counts are distributed for each sample
colors <- ifelse(
  grepl("^BR", names(countData)), "royalblue1",
  ifelse(
    grepl("^MAN", names(countData)), "#35B779FF",
    ifelse(
      grepl("^SB", names(countData)), "#F8766D",
      ifelse(grepl("^TEM", names(countData)), "#440154FF", "black"))))

boxplot(log2(countData), col = colors, ylab = "log2(gene expression)", notch = F, varwidth = F, cex = 0.2, las = 2) 

### Create DeSeq object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData), colData = sampleInfo, design = ~ condition) ## Define design 

### Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

### Get normalized counts
dds <- estimateSizeFactors(dds)

normCounts <- counts(dds, normalized = TRUE) %>% data.frame()
normCounts$target_id <- rownames(normCounts)
normCounts <- normCounts[, c("target_id", names(normCounts)[-ncol(normCounts)])]
rownames(normCounts) <- NULL
write.table(normCounts, file = "./results/normalized.counts.txt", sep = "\t", quote = F, col.names = NA) 

### Principal components analysis (PCA)
### Transform data
vsd <- vst(dds, blind = TRUE) 

### Input is matrix of transformed values
vsdMatrix <- assay(vsd)

pca <- prcomp(t(vsdMatrix), scale. = T)
percentVar <- pca$sdev^2 / sum(pca$sdev^2)

### Create data frame with metadata and PC values for input to ggplot
df <- cbind(sampleInfo, pca$x)
df$site <- factor(df$site, levels = c("Temae", "Beach Rock", "Mangrove", "Sand Beach"))
df$day <- factor(df$day, levels = c("D0", "D1", "D3", "D8", "Juvenile"))
#df$condition <- factor(df$condition, levels = c("Temae_D0", "Beach Rock_D1", "Beach Rock_D3", "Beach Rock_Juvenile",
#                                                "Mangrove_D1", "Mangrove_D3", "Mangrove_D8", "Mangrove_Juvenile",
#                                                "Sand Beach_D1", "Sand Beach_D3", "Sand Beach_D8", "Sand Beach_Juvenile"))

### Choose color/shape scale accordingly
colorsSite <- c("Temae" = "#440154FF", "Beach Rock" = "royalblue1", "Mangrove" = "#35B779FF", "Sand Beach" = "#F8766D")
shapesSite <- c("Temae" = 16, "Beach Rock" = 17, "Mangrove" = 15, "Sand Beach" = 5)
colorsDay <- c("D0" = "#91D1C2FF", "D1" = "#F39B7FFF", "D3" = "#FDE725FF", "D8" = "#9632B8FF", "Juvenile" = "darkgrey")
shapesDay <- c("D0" = 16, "D1" = 17, "D3" = 15, "D8" = 5, "Juvenile" = 6)

### Plot PC1 vs PC2
ggplot(df, aes(x = PC1, y = PC2, color = site, group = site, shape = day)) + 
  geom_point(size = 4) + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "lightgrey") + 
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color = "lightgrey")+
  scale_color_manual(values = colorsSite) + 
  scale_shape_manual(values = shapesDay) +
  #geom_text_repel(aes(label = rownames(df)), size = 2.5) + 
  stat_ellipse(lwd = 0.2) +
  xlab(sprintf("PC1 (%.2f%%)", percentVar[1] * 100)) +
  ylab(sprintf("PC2 (%.2f%%)", percentVar[2] * 100)) +
  theme_bw() +
  theme(axis.title = element_text(colour = "black", size = 13),
        axis.text = element_text(colour = "black", size = 12),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.8)) + 
  guides(color = guide_legend(title = "Site"), shape = guide_legend(title = "Day"))

### Plot PC2 vs PC3
ggplot(df, aes(x = PC2, y = PC3, color = site, group = day, shape = day)) + 
  geom_point(size = 4) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.5, color = "lightgrey") + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "lightgrey")+
  scale_color_manual(values = colorsSite) + 
  scale_shape_manual(values = shapesDay) +
  #geom_text_repel(aes(label = rownames(df)), size = 2.5) + 
  stat_ellipse(lwd = 0.2) +
  xlab(sprintf("PC2 (%.2f%%)", percentVar[2] * 100)) +
  ylab(sprintf("PC3 (%.2f%%)", percentVar[3] * 100)) +
  theme_bw() +
  theme(axis.title = element_text(colour = "black", size = 13),
        axis.text = element_text(colour = "black", size = 12),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.8)) +
  guides(color = guide_legend(title = "Site"), shape = guide_legend(title = "Day"))

#Do PERMANOVA to see differences
install.packages("vegan")
library(vegan)

install.packages("RVAideMemoire")
library("RVAideMemoire")

dis <- vegdist(t(vsdMatrix), method = "euclidean", na.rm = TRUE) 

day <- factor(c(rep(1,5), rep(2,15), rep(3,13), rep(4,10), rep(5,7)), labels = c("D0", "D1", "D3", "D8", "juvenile"))

site <- factor(c(rep(1,4), rep(2,3)), labels = c("Beach Rock", "Mangrove"))

adonis2(dis ~ site, data = D1_sampleInfo, permutations = 10000)

pairwise.perm.manova(dis, site, p.method = "fdr", nperm = 10000)

### Hierarchical Clustering
### Compute pairwise correlation values
vsdCor <- cor(vsdMatrix) #first compute pairwise correlation values

### Plot correlation values as heatmap
sampleInfo$site <- factor(sampleInfo$site, levels = unique(sampleInfo$site))

pheatmap(vsdCor, 
         annotation = sampleInfo[, c("site", "day")],
         scale = "row", 
         clustering_method = "complete",
         cluster_rows = T, cluster_cols = F,
         col = colorRampPalette(c("grey40", "white", "blue"))(256), 
         border_color = NA, 
         fontsize = 6, 
         angle_col = 45, 
         annotation_colors = list(site = colorsSite[unique(sampleInfo$site)],
                                  day = colorsDay[unique(sampleInfo$day)]))

### Compute Euclidean distance between samples
dists	= dist(t(assay(vsd)))

### Perform clustering with hclust
hc <- hclust(dists, method = "complete") ## you can also do ward.D2
dhc <- as.dendrogram(hc)

### Plot dendrogram
### Create function to add attributes to each leaf
colLab <- function(n) {
  if (is.leaf(n)) {
    label <- attributes(n)$label
    day <- sampleInfo$day[match(label, rownames(sampleInfo))]
    day_color <- ifelse(day == "D0", "#91D1C2FF", 
                        ifelse(day == "D1", "#F39B7FFF", 
                               ifelse(day == "D3", "#FDE725FF", 
                                      ifelse(day == "D8", "#9632B8FF", 
                                             ifelse(day == "Juvenile", "darkgrey", NA)))))
                        attr(n, "nodePar") <- list(cex = 1.5, lab.cex = 1, pch = 20, col = day_color, lab.col = day_color, lab.font = 1)
                        }
  return(n)
}

### Apply function to dendrogram
dL <- dendrapply(dhc, colLab)

### Plot the dendrogram with the specified attributes
plot(dL, cex = 0.7, main = "Structure of the Population")

### Create legend 
legend("topright", 
       legend = c("D0", "D1", "D3", "D8", "Juvenile"),
       col = c("#91D1C2FF", "#F39B7FFF", "#FDE725FF", "#9632B8FF", "darkgrey"),
       pch = 20,
       bty = "n",
       pt.cex = 1.5,
       cex = 0.8,
       text.col = "black",
       horiz = FALSE,
       inset = c(0, 0.1))

### Plot heatmap with top 1000 genes with highest variance
### Select top 1000 genes with highest variance

vsdMatrix_D1D3D8 <- vsdMatrix[, c(6:13, 18:32, 36:50)] #D1, D3 & D8
vsdMatrix_D1 <- vsdMatrix[, c(1:10, 18:22, 36:40)] #D0 & D1
vsdMatrix_D3 <- vsdMatrix[, c(11:13, 23:27, 41:45)] #D3
vsdMatrix_D1D3 <- vsdMatrix[, c(6:13, 18:27, 36:45)] #D1 & D3
vsdMatrix_D8 <- vsdMatrix[, c(28:32, 46:50)] #D8
vsdMatrix_D3D8 <- vsdMatrix[, c(11:13, 23:32, 41:50)] #D3 & D8
vsdMatrix_JUV <- vsdMatrix[, c(14:17, 33:35)] #Juveniles
vsdMatrix_D8JUV <- vsdMatrix[, c(14:17, 28:35, 46:50)] #D8 & Juveniles
vsdMatrix_MAN <- vsdMatrix[, c(18:35)] #Mangrove

topVarGenes <- vsdMatrix_D1D3D8[head(order(-rowVars(vsdMatrix_red)), 1000), ] ## using all genes
topVarGenes <- vsdMatrix_JUV[head(order(-rowVars(vsdMatrix_JUV)), 1000), ] ## specify subset of samples you want to check

### Plot heatmap using Pearson's correlation for the dendrogram and using denextend for colors
Rowv_Pear  <- hclust(as.dist(1-cor(t(topVarGenes)))) %>% as.dendrogram %>% set("branches_k_color", k = 12) %>% set("branches_lwd", 1) %>% ladderize ## Choose how many colored branches you want to show
Colv  <- topVarGenes %>% t %>% dist %>% hclust %>% as.dendrogram %>% set("branches_k_color", k = 2) %>% set("branches_lwd", 2) %>% ladderize

### Plot heatmap with cluster labels
heatmap <- heatmap.2(topVarGenes,
                     col = colorRampPalette(c("blue", "white", "red"))(256), 
                     Rowv = Rowv_Pear, Colv = Colv, 
                     scale = "row", dendrogram = "both", trace = "none", density.info = "none", 
                     srtCol = 45, lwid = c(1, 5), lhei = c(2, 7), 
                     key = TRUE, keysize = 1.3)

### Get genes in same order as shown in heatmap
heatmapLabels <- labels(Rowv_Pear) %>% as.data.frame() 
colnames(heatmapLabels)[colnames(heatmapLabels) == "."] <- "target_id"

### Get clusters from heatmap
rowClusters <- as.data.frame(cutree(Rowv_Pear, k = 12)) %>% rownames_to_column("target_id") 
colnames(rowClusters)[colnames(rowClusters) == "cutree(Rowv_Pear, k = 12)"] <- "cluster"

### Assign cluster information to each gene
geneClusterInfo <- heatmapLabels %>%
  left_join(rowClusters, by = "target_id") 
write.table(geneClusterInfo, file = "./GenesHeatmapTop1000VarGenesPearson_Juveniles.txt", sep = "\t", quote = F, col.names = NA) 

### Verify how many genes there are for each cluster
geneClusterInfo %>%
  group_by(cluster) %>%
  summarise(count = n())

######################################2. Differential Expression Analysis######################################
### See https://github.com/hbctraining/DGE_workshop_salmon/blob/master/lessons/01_DGE_setup_and_overview.md for reference
### Run analysis
dds <- DESeq(dds)

############################2.1 Differences across environments############################
##############2.1.1 Sand Beach vs Beach Rock ##############
### Create contrasts 
### contrast <- c("condition", "level_to_compare", "base_level")
contrastD1 <- c("condition", "Sand Beach_D1", "Beach Rock_D1")
contrastD3 <- c("condition", "Sand Beach_D3", "Beach Rock_D3")

### Build results table
resTableSBD1vsBRD1 <- results(dds, contrast = contrastD1, alpha = 0.05) 
resTableSBD3vsBRD3 <- results(dds, contrast = contrastD3, alpha = 0.05) 

### Apply fold change shrinkage
resTableSBD1vsBRD1 <- lfcShrink(dds, contrast = contrastD1, res = resTableSBD1vsBRD1, type = "ashr") 
resTableSBD3vsBRD3 <- lfcShrink(dds, contrast = contrastD3, res = resTableSBD3vsBRD3, type = "ashr") 

### Extract significant DEGs
### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- log2(1.5)

resTable_tbSBD1vsBRD1 <- resTableSBD1vsBRD1 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbSBD1vsBRD1, file = "./SBD1vsBRD1.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

resTable_tbSBD3vsBRD3 <- resTableSBD3vsBRD3 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbSBD3vsBRD3, file = "./SBD3vsBRD3.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve DEGs
sigSBD1vsBRD1 <- resTable_tbSBD1vsBRD1 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigSBD1vsBRD1, file = "./SBD1vsBRD1.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

sigSBD3vsBRD3 <- resTable_tbSBD3vsBRD3 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigSBD3vsBRD3, file = "./SBD3vsBRD3.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve up-regulated genes
upSBD1vsBRD1 <- subset(sigSBD1vsBRD1, log2FoldChange >= lfc.cutoff)
write.table(upSBD1vsBRD1, file = "./SBD1vsBRD1.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

upSBD3vsBRD3 <- subset(sigSBD3vsBRD3, log2FoldChange >= lfc.cutoff)
write.table(upSBD3vsBRD3, file = "./SBD3vsBRD3.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve down-regulated genes
downSBD1vsBRD1 <- subset(sigSBD1vsBRD1, log2FoldChange < lfc.cutoff)
write.table(downSBD1vsBRD1, file = "./SBD1vsBRD1.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

downSBD3vsBRD3 <- subset(sigSBD3vsBRD3, log2FoldChange < lfc.cutoff)
write.table(downSBD3vsBRD3, file = "./SBD3vsBRD3.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

#resTable_tbSBD1vsBRD1 <- resTable_tbSBD1vsBRD1 %>% 
#  mutate(threshold_SBD1vsBRD1 = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
#resTable_tbSBD3vsBRD3 <- resTable_tbSBD3vsBRD3 %>% 
#  mutate(threshold_SBD3vsBRD3 = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

##############2.1.2 Mangrove vs Beach Rock ##############
### Create contrasts 
contrastD1 <- c("condition", "Mangrove_D1", "Beach Rock_D1")
contrastD3 <- c("condition", "Mangrove_D3", "Beach Rock_D3")
contrastJUV <- c("condition", "Mangrove_Juvenile", "Beach Rock_Juvenile")

### Build results table
resTableMAND1vsBRD1 <- results(dds, contrast = contrastD1, alpha = 0.05) 
resTableMAND3vsBRD3 <- results(dds, contrast = contrastD3, alpha = 0.05) 
resTableMANJUVvsBRJUV <- results(dds, contrast = contrastJUV, alpha = 0.05) 

### Apply fold change shrinkage
resTableMAND1vsBRD1 <- lfcShrink(dds, contrast = contrastD1, res = resTableMAND1vsBRD1, type = "ashr") 
resTableMAND3vsBRD3 <- lfcShrink(dds, contrast = contrastD3, res = resTableMAND3vsBRD3, type = "ashr") 
resTableMANJUVvsBRJUV <- lfcShrink(dds, contrast = contrastJUV, res = resTableMANJUVvsBRJUV, type = "ashr") 

### Extract significant DEGs
resTable_tbMAND1vsBRD1 <- resTableMAND1vsBRD1 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbMAND1vsBRD1, file = "./MAND1vsBRD1.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

resTable_tbMAND3vsBRD3 <- resTableMAND3vsBRD3 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbMAND3vsBRD3, file = "./MAND3vsBRD3.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

resTable_tbMANJUVvsBRJUV <- resTableMANJUVvsBRJUV %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbMANJUVvsBRJUV, file = "./MANJUVvsBRJUV.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve DEGs
sigMAND1vsBRD1 <- resTable_tbMAND1vsBRD1 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigMAND1vsBRD1, file = "./MAND1vsBRD1.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

sigMAND3vsBRD3 <- resTable_tbMAND3vsBRD3 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigMAND3vsBRD3, file = "./MAND3vsBRD3.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

sigMANJUVvsBRJUV <- resTable_tbMANJUVvsBRJUV %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigMANJUVvsBRJUV, file = "./MANJUVvsBRJUV.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve up-regulated genes
upMAND1vsBRD1 <- subset(sigMAND1vsBRD1, log2FoldChange >= lfc.cutoff)
write.table(upMAND1vsBRD1, file = "./MAND1vsBRD1.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

upMAND3vsBRD3 <- subset(sigMAND3vsBRD3, log2FoldChange >= lfc.cutoff)
write.table(upMAND3vsBRD3, file = "./MAND3vsBRD3.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

upMANJUVvsBRJUV <- subset(sigMANJUVvsBRJUV, log2FoldChange >= lfc.cutoff)
write.table(upMANJUVvsBRJUV, file = "./MANJUVvsBRJUV.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve down-regulated genes
downMAND1vsBRD1 <- subset(sigMAND1vsBRD1, log2FoldChange < lfc.cutoff)
write.table(downMAND1vsBRD1, file = "./MAND1vsBRD1.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

downMAND3vsBRD3 <- subset(sigMAND3vsBRD3, log2FoldChange < lfc.cutoff)
write.table(downMAND3vsBRD3, file = "./MAND3vsBRD3.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

downMANJUVvsBRJUV <- subset(sigMANJUVvsBRJUV, log2FoldChange < lfc.cutoff)
write.table(downMANJUVvsBRJUV, file = "./MANJUVvsBRJUV.down.txt", sep = "\t", quote = FALSE, col.names = NA) 


#resTable_tbMAND1vsBRD1 <- resTable_tbMAND1vsBRD1 %>% 
#    mutate(threshold_MAND1vsBRD1 = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
#resTable_tbMAND3vsBRD3 <- resTable_tbMAND3vsBRD3 %>% 
#  mutate(threshold_MAND3vsBRD3 = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
#resTable_tbMANJUVvsBRJUV <- resTable_tbMANJUVvsBRJUV %>% 
#  mutate(threshold_MANJUVvsBRJUV = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
  
##############2.1.3 Mangrove vs Sand Beach ##############
### Create contrasts 
contrastD1 <- c("condition", "Mangrove_D1", "Sand Beach_D1")
contrastD3 <- c("condition", "Mangrove_D3", "Sand Beach_D3")
contrastD8 <- c("condition", "Mangrove_D8", "Sand Beach_D8")

### Build results table
resTableMAND1vsSBD1 <- results(dds, contrast = contrastD1, alpha = 0.05) 
resTableMAND3vsSBD3 <- results(dds, contrast = contrastD3, alpha = 0.05) 
resTableMAND8vsSBD8 <- results(dds, contrast = contrastD8, alpha = 0.05) 

### Apply fold change shrinkage
resTableMAND1vsSBD1 <- lfcShrink(dds, contrast = contrastD1, res = resTableMAND1vsSBD1, type = "ashr") 
resTableMAND3vsSBD3 <- lfcShrink(dds, contrast = contrastD3, res = resTableMAND3vsSBD3, type = "ashr") 
resTableMAND8vsSBD8 <- lfcShrink(dds, contrast = contrastD8, res = resTableMAND8vsSBD8, type = "ashr") 

### Extract significant DEGs
resTable_tbMAND1vsSBD1 <- resTableMAND1vsSBD1 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbMAND1vsSBD1, file = "./MAND1vsSBD1.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

resTable_tbMAND3vsSBD3 <- resTableMAND3vsSBD3 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbMAND3vsSBD3, file = "./MAND3vsSBD3.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

resTable_tbMAND8vsSBD8 <- resTableMAND8vsSBD8 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbMAND8vsSBD8, file = "./MAND8vsSBD8.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve DEGs
sigMAND1vsSBD1 <- resTable_tbMAND1vsSBD1 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigMAND1vsSBD1, file = "./MAND1vsSBD1.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

sigMAND3vsSBD3 <- resTable_tbMAND3vsSBD3 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigMAND3vsSBD3, file = "./MAND3vsSBD3.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

sigMAND8vsSBD8 <- resTable_tbMAND8vsSBD8 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigMAND8vsSBD8, file = "./MAND8vsSBD8.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve up-regulated genes
upMAND1vsSBD1 <- subset(sigMAND1vsSBD1, log2FoldChange >= lfc.cutoff)
write.table(upMAND1vsSBD1, file = "./MAND1vsSBD1.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

upMAND3vsSBD3 <- subset(sigMAND3vsSBD3, log2FoldChange >= lfc.cutoff)
write.table(upMAND3vsSBD3, file = "./MAND3vsSBD3.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

upMAND8vsSBD8 <- subset(sigMAND8vsSBD8, log2FoldChange >= lfc.cutoff)
write.table(upMAND8vsSBD8, file = "./MAND8vsSBD8.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve down-regulated genes
downMAND1vsSBD1 <- subset(sigMAND1vsSBD1, log2FoldChange < lfc.cutoff)
write.table(downMAND1vsSBD1, file = "./MAND1vsSBD1.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

downMAND3vsSBD3 <- subset(sigMAND3vsSBD3, log2FoldChange < lfc.cutoff)
write.table(downMAND3vsSBD3, file = "./MAND3vsSBD3.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

downMAND8vsSBD8 <- subset(sigMAND8vsSBD8, log2FoldChange < lfc.cutoff)
write.table(downMAND8vsSBD8, file = "./MAND8vsSBD8.down.txt", sep = "\t", quote = FALSE, col.names = NA) 


#resTable_tbMAND1vsSBD1 <- resTable_tbMAND1vsSBD1 %>% 
#    mutate(threshold_MAND1vsSBD1 = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
#resTable_tbMAND3vsSBD3 <- resTable_tbMAND3vsSBD3 %>% 
#  mutate(threshold_MAND3vsSBD3 = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
#resTable_tbMAND8vsSBD8 <- resTable_tbMAND8vsSBD8 %>% 
#  mutate(threshold_MAND8vsSBD8 = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

############################2.2 Differences across time points############################
##############2.2.1 D1 vs D0 ##############
### Create contrasts 
contrastBR <- c("condition", "Beach Rock_D1", "Temae_D0")
contrastSB <- c("condition", "Sand Beach_D1", "Temae_D0")
contrastMAN <- c("condition", "Mangrove_D1", "Temae_D0")

### Build results table
resTableBRD1vsTEMD0 <- results(dds, contrast = contrastBR, alpha = 0.05) 
resTableSBD1vsTEMD0 <- results(dds, contrast = contrastSB, alpha = 0.05) 
resTableMAND1vsTEMD0 <- results(dds, contrast = contrastMAN, alpha = 0.05) 

### Apply fold change shrinkage
resTableBRD1vsTEMD0 <- lfcShrink(dds, contrast = contrastBR, res = resTableBRD1vsTEMD0, type = "ashr") 
resTableSBD1vsTEMD0 <- lfcShrink(dds, contrast = contrastSD, res = resTableSBD1vsTEMD0, type = "ashr") 
resTableMAND1vsTEMD0 <- lfcShrink(dds, contrast = contrastMAN, res = resTableMAND1vsTEMD0, type = "ashr") 

### Extract significant DEGs
resTable_tbBRD1vsTEMD0 <- resTableBRD1vsTEMD0 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbBRD1vsTEMD0, file = "./BRD1vsTEMD0.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

resTable_tbSBD1vsTEMD0 <- resTableSBD1vsTEMD0 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbSBD1vsTEMD0, file = "./SBD1vsTEMD0.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

resTable_tbMAND1vsTEMD0 <- resTableMAND1vsTEMD0 %>% 
  data.frame() %>% 
  rownames_to_column(var = "target_id") ## This is table with all genes
write.table(resTable_tbMAND1vsTEMD0, file = "./MAND1vsTEMD0.allGenes.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve DEGs
sigBRD1vsTEMD0 <- resTable_tbBRD1vsTEMD0 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigBRD1vsTEMD0, file = "./BRD1vsTEMD0.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

sigSBD1vsTEMD0 <- resTable_tbSBD1vsTEMD0 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigSBD1vsTEMD0, file = "./SBD1vsTEMD0.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

sigMAND1vsTEMD0 <- resTable_tbMAND1vsTEMD0 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
write.table(sigMAND1vsTEMD0, file = "./MAND1vsTEMD0.DEGs.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve up-regulated genes
upBRD1vsTEMD0 <- subset(sigBRD1vsTEMD0, log2FoldChange >= lfc.cutoff)
write.table(upBRD1vsTEMD0, file = "./BRD1vsTEMD0.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

upSBD1vsTEMD0 <- subset(sigSBD1vsTEMD0, log2FoldChange >= lfc.cutoff)
write.table(upSBD1vsTEMD0, file = "./SBD1vsTEMD0.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

upMAND1vsTEMD0 <- subset(sigMAND1vsTEMD0, log2FoldChange >= lfc.cutoff)
write.table(upMAND1vsTEMD0, file = "./MAND1vsTEMD0.up.txt", sep = "\t", quote = FALSE, col.names = NA) 

### Retrieve down-regulated genes
downBRD1vsTEMD0 <- subset(sigBRD1vsTEMD0, log2FoldChange < lfc.cutoff)
write.table(downBRD1vsTEMD0, file = "./BRD1vsTEMD0.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

downSBD1vsTEMD0 <- subset(sigSBD1vsTEMD0, log2FoldChange < lfc.cutoff)
write.table(downSBD1vsTEMD0, file = "./SBD1vsTEMD0.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

downMAND1vsTEMD0 <- subset(sigMAND1vsTEMD0, log2FoldChange < lfc.cutoff)
write.table(downMAND1vsTEMD0, file = "./MAND1vsTEMD0.down.txt", sep = "\t", quote = FALSE, col.names = NA) 

############################2.3 Plot all DEGs############################
### Make file (in Excel) with DEGs for all comparisons 
results <- read.table("./DEGs.comparisons.txt", header = TRUE, sep = "\t")

## Alternative plot
ggplot(results, aes(comparison, y = degs, fill = direction)) + 
  geom_bar(data = subset(results, direction == "upRegulated"), aes(y = degs), stat = "identity", width = 0.8) +
  geom_bar(data = subset(results, direction == "downRegulated"), aes(y = -degs), stat = "identity", width = 0.8) + 
  scale_fill_manual(values = c("upRegulated" = "red", "downRegulated" = "blue"), labels = c("upRegulated" = "up-regulated", "downRegulated" = "down-regulated")) +
  geom_hline(yintercept = 0, colour = "white", size = 1) +
  scale_y_continuous(limits = c(-4500, 3500), breaks = seq(-4500, 3500, by = 1000), labels = abs) +
  labs(x = "", y = "number of DE genes") +
  geom_text(data = subset(results, direction == "upRegulated"), 
            aes(comparison, y = degs, label = abs(degs)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  geom_text(data = subset(results, direction == "downRegulated"), 
            aes(comparison, y = -degs, label = abs(degs)),
            position = position_dodge(width = 0.9), vjust = 1.5,  size = 3) +
  coord_cartesian(ylim = c(-4500, 3500)) +
  theme_bw() +
  theme(axis.title = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        axis.text.x = element_text(colour = "black", size = 8, angle = 45),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.5),
        legend.position = c(0.59, 0.89),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8)) 

######################################3. Plot PCAs and heatmaps for gene sets######################################
### Read files with gene sets
dataFiles <- c(
  ## Pigmentation genes
  #"manini.all.melanophores.txt",
  #"manini.all.xanthophores.txt",
  "manini.components.melanosomes.txt",  
  "manini.iridophores.txt",
  "manini.leucophores.txt",
  "manini.melanocytes.development.txt",
  "manini.melanogenesis.regulation.txt",
  "manini.melanophores.development.txt",
  "manini.melanosome.biogenesis.txt",
  "manini.melanosome.transport.txt", 
  "manini.pigment.cell.differentiation.txt",
  "manini.xanthophore.pteridine.synthesis.txt",
  "manini.xanthophores.development.txt",
  ## Other genes
  #"manini.all.digestion.txt",
  #"manini.all.intestinal.development.txt",
  #"manini.all.stomach.txt",
  "manini.aging.txt",
  "manini.appetite.txt",
  "manini.beta.oxidation.txt",
  "manini.cholesterol.biosynthesis.txt",
  "manini.corticoids.txt",
  "manini.corticosteroids.txt",
  "manini.digestion.glucosidases.txt",
  "manini.digestion.lipases.txt",
  "manini.digestion.proteases.txt",
  "manini.fatty.acid.metabolism.txt",
  "manini.gaba.pathway.txt",
  "manini.glycolysis.txt",
  "manini.intestinal.development.txt",
  "manini.intestinal.development.stefano.txt",
  "manini.ion.channels.txt",
  "manini.krebs.cycle.txt",
  "manini.lactic.fermentation.txt",
  "manini.msh.txt",
  "manini.neurotransmitter.markers.txt",
  "manini.osmoregulation.txt",
  "manini.ossification.txt",
  "manini.pdh.complex.txt",
  "manini.phototransduction.txt",
  "manini.prolactin.txt",
  "manini.retinoic.acid.txt",
  "manini.ROS.txt",
  "manini.stomach.development.txt",
  "manini.stomach.specification.txt",
  "manini.th.pathway.txt",
  "manini.vision.txt"
)

for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("/Users/marleenklann/Desktop/manini/gene_lists/", file), header = TRUE, sep = "\t")
  
  merged <- merge(sigMAND3vsBRD3, data, by = "target_id")
  if("gene_id" %in% names(merged)) {
    # Move gene_id to the first position
    cols <- c("gene_id", setdiff(names(merged), "gene_id"))
    merged <- merged[, cols]
  }
  outputFile <- paste0("./results/sig_MAND3vsBRD3_", gsub("manini\\.", "", file))
  write.table(merged, file = outputFile, sep = "\t", row.names = FALSE, quote = FALSE)
} 

### File should look like this (2 columns):
#gene_id	target_id
#acd	jg4869
#cct5	jg21290
#ctc1	jg30195

### Retrieve data for specific comparisons to plot
##D1
D1 <- countData %>% select(contains("D1"))
D1 <- D1 %>% select(sort(names(D1)))
sampleInfoD1 <- sampleInfo %>% filter(day == "D1")
sampleInfoD1 <- sampleInfoD1 %>% arrange(rownames(sampleInfoD1))

ddsD1 <- DESeqDataSetFromMatrix(countData = as.matrix(D1), colData = sampleInfoD1, design = ~ condition)
keep <- rowSums(counts(ddsD1)) >= 10
ddsD1 <- ddsD1[keep,]
ddsD1 <- estimateSizeFactors(ddsD1)
vsdD1 <- vst(ddsD1, blind = TRUE) 
vsdMatrixD1 <- assay(vsdD1)

ddsD1 <- estimateSizeFactors(ddsD1)
normCountsD1 <- counts(ddsD1, normalized = TRUE) %>% data.frame()
normCountsD1$target_id <- rownames(normCountsD1)
normCountsD1 <- normCountsD1[, c("target_id", names(normCountsD1)[-ncol(normCountsD1)])]
rownames(normCountsD1) <- NULL

##D3
D3 <- countData %>% select(contains("D3"))
D3 <- D3 %>% select(sort(names(D3)))
sampleInfoD3 <- sampleInfo %>% filter(day == "D3")
sampleInfoD3 <- sampleInfoD3 %>% arrange(rownames(sampleInfoD3))

ddsD3 <- DESeqDataSetFromMatrix(countData = as.matrix(D3), colData = sampleInfoD3, design = ~ condition)
keep <- rowSums(counts(ddsD3)) >= 10
ddsD3 <- ddsD3[keep,]
ddsD3 <- estimateSizeFactors(ddsD3)
vsdD3 <- vst(ddsD3, blind = TRUE) 
vsdMatrixD3 <- assay(vsdD3)

ddsD3 <- estimateSizeFactors(ddsD3)
normCountsD3 <- counts(ddsD3, normalized = TRUE) %>% data.frame()
normCountsD3$target_id <- rownames(normCountsD3)
normCountsD3 <- normCountsD3[, c("target_id", names(normCountsD3)[-ncol(normCountsD3)])]
rownames(normCountsD3) <- NULL

##D8
D8 <- countData %>% select(contains("D8"))
D8 <- D8 %>% select(sort(names(D8)))
sampleInfoD8 <- sampleInfo %>% filter(day == "D8")
sampleInfoD8 <- sampleInfoD8 %>% arrange(rownames(sampleInfoD8))

ddsD8 <- DESeqDataSetFromMatrix(countData = as.matrix(D8), colData = sampleInfoD8, design = ~ condition)
keep <- rowSums(counts(ddsD8)) >= 10
ddsD8 <- ddsD8[keep,]
ddsD8 <- estimateSizeFactors(ddsD8)
vsdD8 <- vst(ddsD8, blind = TRUE) 
vsdMatrixD8 <- assay(vsdD8)

ddsD8 <- estimateSizeFactors(ddsD8)
normCountsD8 <- counts(ddsD8, normalized = TRUE) %>% data.frame()
normCountsD8$target_id <- rownames(normCountsD8)
normCountsD8 <- normCountsD8[, c("target_id", names(normCountsD8)[-ncol(normCountsD8)])]
rownames(normCountsD8) <- NULL

##Juveniles
JUV <- countData %>% select(contains("JUV"))
JUV <- JUV %>% select(sort(names(JUV)))
sampleInfoJUV <- sampleInfo %>% filter(day == "Juvenile")
sampleInfoJUV <- sampleInfoJUV %>% arrange(rownames(sampleInfoJUV))

ddsJUV <- DESeqDataSetFromMatrix(countData = as.matrix(JUV), colData = sampleInfoJUV, design = ~ condition)
keep <- rowSums(counts(ddsJUV)) >= 10
ddsJUV <- ddsJUV[keep,]
ddsJUV <- estimateSizeFactors(ddsJUV)
vsdJUV <- vst(ddsJUV, blind = TRUE) 
vsdMatrixJUV <- assay(vsdJUV)

ddsJUV <- estimateSizeFactors(ddsJUV)
normCountsJUV <- counts(ddsJUV, normalized = TRUE) %>% data.frame()
normCountsJUV$target_id <- rownames(normCountsJUV)
normCountsJUV <- normCountsJUV[, c("target_id", names(normCountsJUV)[-ncol(normCountsJUV)])]
rownames(normCountsJUV) <- NULL

##D0 and D1
D0vsD1 <- countData %>% select(contains("D0") | contains("D1")) 
D0vsD1 <- D0vsD1 %>% select(sort(names(D0vsD1)))
sampleInfoD0vsD1 <- sampleInfo %>% filter(day == "D0" | day == "D1")
sampleInfoD0vsD1 <- sampleInfoD0vsD1 %>% arrange(rownames(sampleInfoD0vsD1))

ddsD0vsD1 <- DESeqDataSetFromMatrix(countData = as.matrix(D0vsD1), colData = sampleInfoD0vsD1, design = ~ condition)
keep <- rowSums(counts(ddsD0vsD1)) >= 10
ddsD0vsD1 <- ddsD0vsD1[keep,]
ddsD0vsD1 <- estimateSizeFactors(ddsD0vsD1)
vsdD0vsD1 <- vst(ddsD0vsD1, blind = TRUE) 
vsdMatrixD0vsD1 <- assay(vsdD0vsD1)

ddsD0vsD1 <- estimateSizeFactors(ddsD0vsD1)
normCountsD0vsD1 <- counts(ddsD0vsD1, normalized = TRUE) %>% data.frame()
normCountsD0vsD1$target_id <- rownames(normCountsD0vsD1)
normCountsD0vsD1 <- normCountsD0vsD1[, c("target_id", names(normCountsD0vsD1)[-ncol(normCountsD0vsD1)])]
rownames(normCountsD0vsD1) <- NULL

### Loop through files
for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("../../../gene_lists/", file), header = TRUE, sep = "\t")
  
  ### Filter vsd_matrix
  vsdDf <- as.data.frame(vsdMatrix) ## Make sure to use correct dataset
  vsdDf <- rownames_to_column(vsdDf, var = "target_id") 
  geneCounts <- vsdDf %>% filter(target_id %in% data$target_id)
  row.names(geneCounts) <- geneCounts[, 1]
  geneCounts <- geneCounts[, -1]
  #geneCounts <- geneCounts[, c(1:10, 18:22, 36:40)] #D0 & D1
  #geneCounts <- geneCounts[, c(6:13, 18:27, 36:45)] #D1 & D3
  #geneCounts <- geneCounts[, c(11:13, 23:27, 41:45)] #D3
  #geneCounts <- geneCounts[, c(11:13, 23:32, 41:50)] #D3 & D8
  #geneCounts <- geneCounts[, c(28:32, 46:50)] #D8
  #geneCounts <- geneCounts[, c(14:17, 28:35, 46:50)] #D8 & Juveniles
  #geneCounts <- geneCounts[, c(14:17, 33:35)] #Juveniles
  geneCounts <- geneCounts[, c(6:13, 18:32, 36:50)] #D1 & D3 & D8
  #geneCounts <- geneCounts[, c(1:13, 18:32, 36:50)] #D0 & D1 & D3 & D8
  
  ### Perform PCA clustering based on specific gene set
  pca <- prcomp(t(geneCounts), scale. = T)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  #df <- cbind(sampleInfo, pca$x)
  df <- cbind((sampleInfo %>% filter(!grepl("D0|Juvenile", day))), pca$x) 
  df$site <- factor(df$site, levels = c("Beach Rock", "Mangrove", "Sand Beach"))
  df$day <- factor(df$day, levels = c("D1", "D3", "D8"))
  df$condition <- factor(df$condition, levels = c("Beach Rock_D1", "Beach Rock_D3", #"Beach Rock_Juvenile",
                                                  "Mangrove_D1", "Mangrove_D3", "Mangrove_D8",
                                                  "Sand Beach_D1", "Sand Beach_D3", "Sand Beach_D8"))
  #df <- cbind((sampleInfo %>% filter(!grepl("D0|D1|D3|D8", day))), pca$x) ## Make sure to use correct metadata file
  #df$day <- factor(df$day, levels = c("Juvenile"))
  #df$site <- factor(df$site, levels = c("Beach Rock", "Mangrove"))
  
  ### Plot PC1 vs PC2
  colorsSite <- c("Beach Rock" = "royalblue1", "Mangrove" = "#35B779FF", "Sand Beach" = "#F8766D") #"Temae" = "#440154FF", 
  #shapesSite <- c("Temae" = 16, "Beach Rock" = 17, "Mangrove" = 15, "Sand Beach" = 5) 
  #colorsDay <- c("D1" = "#F39B7FFF", "D3" = "#FDE725FF", "D8" = "#9632B8FF", "D0" = "#91D1C2FF")
  shapesDay <- c("D1" = 17, "D3" = 15, "D8" = 5) #"Juvenile" = 6 "D0" = 16, 
  #shapesDay <- c("Juvenile" = 6)
  
  ### Adjust settings according to data set you want to plot
  pca_plot <- ggplot(df, aes(x = PC1, y = PC2, color = site, group = condition, shape = day)) + 
    geom_point(size = 4) + 
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "lightgrey") + 
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color = "lightgrey")+
    scale_color_manual(values = colorsSite, labels = c("Beach Rock", "Mangrove", "Sand Beach")) + 
    scale_shape_manual(values = shapesDay, labels = c("D1", "D3", "D8")) +
    stat_ellipse(lwd = 0.2) +
    xlab(sprintf("PC1 (%.2f%%)", percentVar[1] * 100)) +
    ylab(sprintf("PC2 (%.2f%%)", percentVar[2] * 100)) +
    theme_bw() +
    theme(
      axis.title = element_text(colour = "black", size = 13),
      axis.text = element_text(colour = "black", size = 12),
      axis.ticks = element_line(colour = "black", size = 0.5),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", size = 0.8)) +
  guides(color = guide_legend(title = "Site"), shape = guide_legend(title = "Day"))

  ### Save each PCA plot
  ggsave(pca_plot, filename = paste0("./figures/PCA_", gsub("manini\\.|\\.txt", "", file), "_D1D3D8_condition", ".pdf"), width = 8, height = 6)
} 

for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("/Users/marleenklann/Desktop/manini/gene_lists/", file), header = TRUE, sep = "\t")
  
  ### Prepare input data
  vsdDf <- as.data.frame(vsdMatrix)
  vsdDf <- rownames_to_column(vsdDf, var = "target_id") 
  geneCounts <- vsdDf %>% filter(target_id %in% data$target_id)
  geneCounts <- merge(geneCounts, data, by.x = "target_id", by.y = "target_id", all.x = TRUE)
  row.names(geneCounts) <- geneCounts$gene_id
  geneCounts <- geneCounts[, c("gene_id", names(vsdDf)[-1])]
  geneCounts <- geneCounts[, -1]
  #geneCounts <- geneCounts[, c(6:13, 18:32, 36:50)] #D1 & D3 & D8
  #column_order <- c("BR", "SB", "MAN")
  #time_points <- c("D1", "D3", "D8")
  #geneCounts <- geneCounts[, order(gsub(".*_(.*?)_.*", "\\1", colnames(geneCounts)), 
  #                                 match(gsub(".*_.*_(.*?)", "\\1", colnames(geneCounts)), time_points))]
  geneCounts <- geneCounts[, c(14:17, 33:35)] #Juveniles
  
  ### Plot heatmap 
  heatmap <- pheatmap(geneCounts, 
                      scale = "row", 
                      col = colorRampPalette(c("blue", "white", "red"))(256),
                      border_color = NA, 
                      #cluster_rows = FALSE, 
                      cluster_cols = FALSE,
                      fontsize_row = 8, angle_col = 45,
                      cellwidth = 20, cellheight = 8.5)
  
  ### Save each plot plot
  ggsave(heatmap, filename = paste0("./figures/heatmap_", gsub("manini\\.|\\.txt", "", file), "_Juveniles", ".pdf"), width = 15, height = 10)
} 

######################################4. Plot Expression Levels######################################
### Run ANOVAs to see differences in expression levels (use normalized counts)

### Loop through all gene sets
for (file in dataFiles) {
  ### Read data
  data <- read.table(file.path("../../../gene_lists/", file), header = TRUE, sep = "\t")
  
  ### Format data frames
  geneCounts <- normCountsD3 %>% inner_join(data, by = "target_id") %>% ## Make sure use correct file
    select(gene_id, everything(), -target_id) %>%
    arrange(gene_id)
  
  geneCounts <- geneCounts %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "Value") %>%
    mutate(site = case_when(grepl("^BR", sample) ~ "Beach Rock", ## Make sure information is correct
                            grepl("^MAN", sample) ~ "Mangrove",
                            grepl("^SB", sample) ~ "Sand Beach",
                            #grepl("^TEM", sample) ~ "Temae D0",
                            TRUE ~ "Other")) %>%
    select(sample, site, gene_id, Value) %>%
    pivot_wider(names_from = gene_id, values_from = Value)

  ### Run ANOVA and save results to a file
  genes <- colnames(geneCounts)[!(colnames(geneCounts) %in% c("sample", "site"))]
  outputFile <- paste0("./anova_expression.levels/anova_D3_", gsub("manini\\.", "", file))
  
  sink(outputFile)
  
  for (gene in genes) {
    formula <- as.formula(paste(gene, "~ site"))
    cat("Gene:", gene, "\n")
    anovaResult <- aov(formula, data = geneCounts)
    anovaSummary <- summary(anovaResult)
    cat("ANOVA Summary:\n", file = outputFile, append = TRUE)
    print(anovaSummary, quote = FALSE, file = outputFile, append = TRUE)
    tukeyResults <- TukeyHSD(anovaResult, "site")
    cat("Tukey HSD Results:\n", file = outputFile, append = TRUE)
    print(tukeyResults, quote = FALSE, file = outputFile, append = TRUE)
    cat("\n", file = outputFile, append = TRUE)
  }
  sink()
}

### Make line plots for time series
## Create empty list
geneCountsList <- list()

## Loop through data files
for (file in dataFiles) {
  ### Read the data
  data <- read.table(file.path("../../../gene_lists/", file), header = TRUE, sep = "\t")
  
  ### Format data frames
  geneCounts <- normCounts[, c(1, 7:14, 19:33, 37:51)] %>% inner_join(data, by = "target_id") %>% ## Make sure to select columns of interest D1, D3, D8
   #geneCounts <- normCounts[, c(1, 7:14, 19:28, 37:46)] %>% inner_join(data, by = "target_id") %>%
  #geneCounts <- normCounts[, c(1, 7:51)] %>% inner_join(data, by = "target_id") %>% ## all samples
    dplyr::select(gene_id, everything(), -target_id) %>%
    arrange(gene_id) %>%
    mutate_at(vars(-gene_id), ~ log2(. + 1))
  
  geneCounts <- geneCounts %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "Value") %>%
    mutate(day = case_when(grepl("D1", sample) ~ "D1",
                           grepl("D3", sample) ~ "D3",
                           grepl("D8", sample) ~ "D8",
                           TRUE ~ "Other")) %>%
    mutate(site = case_when(grepl("^BR", sample) ~ "Beach Rock",
                            grepl("^MAN", sample) ~ "Mangrove",
                           grepl("^SB", sample) ~ "Sand Beach",
                           TRUE ~ "Other")) %>%
    dplyr::select(sample, day, site, gene_id, Value) %>%
    pivot_wider(names_from = gene_id, values_from = Value)
  
  geneCounts$day <- factor(geneCounts$day, levels = c("D1", "D3", "D8"))
  
  ### Save results in the list
  geneCountsList[[file]] <- geneCounts
}

### Generate and save plots for each gene set
for (file in dataFiles) {
  geneCounts <- geneCountsList[[file]]

  ### Extract relevant columns for plotting
  cols_to_plot <- colnames(geneCounts)[!(colnames(geneCounts) %in% c("sample", "day", "site"))]
 
  ### Split the columns into groups of 12 (or depending in how many plots you want)
  col_groups <- split(cols_to_plot, rep(1:(ceiling(length(cols_to_plot) / 12)), each = 12, length.out = length(cols_to_plot)))
  
  ### Generate and save plots for each group 
  for (i in seq_along(col_groups)) {
    cols <- col_groups[[i]]
    
    plots <- list()
    for (col in cols) {
      p <- ggplot(data = geneCounts, aes(x = day, y = !!sym(col), color = site, fill = site, shape = day)) +
        geom_point(position = position_jitter(width = 0.1), size = 1.8, show.legend = FALSE) + # alpha = 0.9, stroke = 0.01
        geom_line(data = geneCounts %>% 
                    group_by(day, site) %>%
                    summarize(mean_value = mean(!!sym(col)), .groups = "drop"),
                  aes(x = day, y = mean_value, group = site, color = site), size = 0.8, show.legend = FALSE) + 
        labs(title = paste(col, ""), x = "", y = "Expression level") +
        scale_shape_manual(values = c("D1" = 17, "D3" = 15, "D8" = 5)) + 
        scale_color_manual(values = c("Beach Rock" = "royalblue1", "Mangrove" = "#35B779FF", "Sand Beach" = "#F8766D")) + 
        theme_classic() + theme(
          plot.title = element_text(size = 9),
          axis.title = element_text(colour = "black", size = 8),
          axis.text = element_text(colour = "black", size = 7),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
          )
      
      plots[[col]] <- p
    }
    
    ### Save plots
    filename <- gsub("manini\\.|\\.txt", "", file)
    pdf(paste0("./figures/expression.lineplots_D1D3D8_", filename, "_", i, ".pdf"))
    grid.arrange(grobs = plots, ncol = 3, nrow = 4)
    dev.off()
  }
}

### Make box plots for JUVENILE comparisons
## Create empty list
geneCountsList <- list()

## Loop through data files
for (file in dataFiles) {
  ### Read the data
  data <- read.table(file.path("../../../gene_lists/", file), header = TRUE, sep = "\t")
  
  ### Format data frames
  geneCounts <- normCounts[, c(1, 15:18, 34:36)]  %>% inner_join(data, by = "target_id") %>%  
    select(gene_id, everything(), -target_id) %>%
    arrange(gene_id) %>%
    mutate_at(vars(-gene_id), ~ log2(. + 1))
  
  geneCounts <- geneCounts %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "Value") %>%
    mutate(day = case_when(grepl("JUV", sample) ~ "JUV",
                           TRUE ~ "Other")) %>%
    mutate(site = case_when(grepl("^BR", sample) ~ "BR",
                            grepl("^MAN", sample) ~ "MAN",
                            TRUE ~ "Other")) %>%
    select(sample, day, site, gene_id, Value) %>%
    pivot_wider(names_from = gene_id, values_from = Value)
  
  geneCounts$day <- factor(geneCounts$day, levels = c("JUV"))
  geneCounts$site <- factor(geneCounts$site, levels = c("BR", "MAN"))
  
  ### Save results in the list
  geneCountsList[[file]] <- geneCounts
}

### Generate and save plots for each gene set
for (file in dataFiles) {
  geneCounts <- geneCountsList[[file]]
  
  ### Extract relevant columns for plotting
  cols_to_plot <- colnames(geneCounts)[!(colnames(geneCounts) %in% c("sample", "day", "site"))]
  
  ### Split the columns into groups of 12 (or depending in how many plots you want)
  col_groups <- split(cols_to_plot, rep(1:(ceiling(length(cols_to_plot) / 12)), each = 12, length.out = length(cols_to_plot)))
  
  ### Generate and save plots for each group 
  for (i in seq_along(col_groups)) {
    cols <- col_groups[[i]]
    
    plots <- list()
    for (col in cols) {
      p <- ggplot(data = geneCounts, aes(x = day, y = !!sym(col), color = site, fill = site)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA, size = 0.4, width = 0.5, show.legend = FALSE) + 
        geom_point(position = position_jitterdodge(), size = 1, show.legend = FALSE) +
        labs(title = paste(col, ""), x = "", y = "Expression level") +
        scale_color_manual(values = c("BR" = "royalblue1", "MAN" = "#35B779FF")) +
        scale_fill_manual(values = c("BR" = "royalblue1", "MAN" = "#35B779FF")) +
        theme_bw() + theme(
          plot.title = element_text(size = 9),
          axis.title = element_text(colour = "black", size = 8),
          axis.text = element_text(colour = "black", size = 7),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8))
      
      plots[[col]] <- p
    }
    
    ### Save plots
    filename <- gsub("manini\\.|\\.txt", "", file)
    pdf(paste0("./figures/expression.boxplots.Juveniles_", filename, "_", i, ".pdf"))
    grid.arrange(grobs = plots, ncol = 3, nrow = 4)
    dev.off()
  }
}
