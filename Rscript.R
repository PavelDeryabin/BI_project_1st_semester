
###########################################################################################################################################

#                                                   RNA-seq analysis for OUABAIN THEME

###########################################################################################################################################

library(Rmisc)
library(xlsx)
library(ggplot2)
library(dplyr)
library(tximeta)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicFeatures)
library(DESeq2)
library(clusterProfiler)
library(DOSE)
library(topGO)
library(enrichplot)
library(pheatmap)
library(genefilter)
library(RColorBrewer)

#### Tximeta quantification 

setwd("...") # setting wd with the folder "quants" and predesigned "coldata.csv" file

dir <- ("quants")
list.files(file.path(dir))

dir2 <- (".")
csvfile <- file.path(dir2, "coldatadf.csv")
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)

str(coldata)
coldata$SampleID <- as.character(coldata$SampleID)
coldata$names <- coldata$SampleGName 
coldata$files <- file.path(dir2, "quants", coldata$names, "quant.sf")
file.exists(coldata$files)
coldata # look at the data

se <- tximeta(coldata) # quantification to transcript-base level
dim(se) # look at the munber of dimentions
head(rownames(se))

gse <- summarizeToGene(se) # summarization to gene-level

keytypes(org.Hs.eg.db) # add nessessary annotations
gse <- addIds(gse, "REFSEQ", gene=TRUE)
gse <- addIds(gse, "ENTREZID", gene=TRUE)
gse <- addIds(gse, "GENENAME", gene=TRUE)
gse <- addIds(gse, "SYMBOL", gene=TRUE)

mcols(gse)

######################## EVA and DE analysis with DESeq2

### Construction of the DESeqDataSet - the DESeq2 custom class object (almost* the same as Bioconductor custom class object SummarizedExperiment)

gse@colData
str(gse@colData$CellType) 
gse@colData$CellType <- as.factor(gse@colData$CellType)
gse@colData$Cells <- as.factor(gse@colData$Cells)
gse@colData$OuabSeno <- as.factor(gse@colData$OuabSeno)
levels(gse@colData$CellType)
levels(gse@colData$Cells)
levels(gse@colData$OuabSeno)
colData(gse)$CellType <- relevel(colData(gse)$CellType, ref = "END-MSCs")
levels(gse@colData$CellType)

levels(gse@colData$Cells) <- c("Young", "Senescent")
levels(gse@colData$CellType) <- c("END-MSCs", "A549", "IMR-90")
levels(gse@colData$OuabSeno) <- c("Sensitive", "Resistant")

names(gse@colData)[names(gse@colData) == 'OuabSeno'] <- 'Senolysis'
levels(gse@colData$Senolysis)

dds <- DESeqDataSet(gse, design = ~ Cells + Senolysis + Cells:Senolysis)
dds # get summary
nrow(dds) # learn total number of rows

### EVA 

dds1 <- dds[ rowSums(counts(dds)) > 5, ] # removing rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds1)
# or
keep <- rowSums(counts(dds) >= 5) >= 4 # minimal filtering to reduce the size of the dataset. We do not need to retain genes 
# if they do not have a count of 5 or more for 4 or more samples 
# as these genes will have no statistical power to detect differences, 
# and no information to compute distances between samples
table(keep)
dds1 <- dds[keep,]

rld <- DESeq2::rlog(dds1, blind = FALSE) # DESeq2 transformation for count data that stabilize the variance across the mean 
class(rld) # DESeqTransform object which is based on the SummarizedExperiment class

DESeq2::plotPCA(rld, intgroup = c("CellType", "Cells"))

# Senescence validation

total <- as.data.frame(dds1@rowRanges)
View(total)
total <- total[,c(6,8)]
total1 <- total[,c(6,10)]
sub1 <- read.xlsx("Cellular senescence GO 0090398.xlsx", sheetIndex = 1)
str(sub1)
str(total)
sub1$ENTREZID <- as.character(sub1$ENTREZID)
merged <- merge(total,sub1,by="ENTREZID")
head(merged)
str(merged)
sub1v <- merged[,2]

rld_sen <- rld[sub1v,]

DESeq2::plotPCA(rld_sen, intgroup = c("CellType", "Cells"))

pcaData <- plotPCA(rld_sen, intgroup = c("CellType", "Cells"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Cells, shape = CellType)) +
  geom_point(size = 2.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(
    plot.title = element_text(size=14, color="black"),
    axis.title.y = element_text(size=14, color="black"),
    axis.title.x = element_text(size=14, color="black"),
    axis.text.y = element_text(size=14),
    axis.text.y.right = element_text(color = "black"),
    axis.text.x = element_text(size=14), 
    legend.text = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 14, colour = "black"),
    legend.position = "right") +
  scale_color_manual(values = c("#ff9a00","#5e5e5e"),  
                     breaks=c("Young","Senescent"),
                     labels=c("Young", "Senescent"))

ggsave(file = "./results_new/PCA.jpeg", width = 4.7, height = 3, dpi = 500)

############ DEG analysis with DESeq2

dds1 <- DESeq(dds1) # running DE analysis (on non-transformed/normilized counts!!!)

resultsNames(dds1) # get terms from the model

res_ouab <- results(dds1, name = "CellsSenescent.SenolysisResistant") # DEGs for interaction

# for now the difference reflects lines vs eMSCs -> *(-1) lfc

mcols(res_ouab, use.names = TRUE) # description of DataFrames
summary(res_ouab) # main statistics

plotMA(res_ouab, ylim = c(-10, 10)) # look at all genes
#abline(h=c(-1,1), col="dodgerblue", lwd=2)

# advanced lcf shrinkage method (more informative visualization and more accurate ranking of genes by effect size - adequate lcf values)
res_ouab_sh <- lfcShrink(dds1, coef = 4, type="apeglm", lfcThreshold=0.667) 

plotMA(res_ouab_sh, ylim = c(-10, 10))
#abline(h=c(-1,1), col="dodgerblue", lwd=2)

### Annotations with different IDs

ens.str <- substr(rownames(res_ouab_sh), 1, 15)
columns(org.Hs.eg.db) # what IDs type are avaliable

res_ouab_sh$symbol <- mapIds(org.Hs.eg.db,
                          keys=ens.str,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

res_ouab_sh$entrez <- mapIds(org.Hs.eg.db,
                          keys=ens.str,
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
res_ouab_sh$genename <- mapIds(org.Hs.eg.db,
                            keys=ens.str,
                            column="GENENAME",
                            keytype="ENSEMBL",
                            multiVals="first")

summary(res_ouab_sh)

############### subsetting 

res_ouab_sh <- res_ouab_sh[order(res_ouab_sh$log2FoldChange, decreasing = T),] # order by lfc
write.csv(as.data.frame(res_ouab_sh), file="./results_new/Supplementary9.csv") # ALL BH-sig DEGs, ordered by lfc

########################################################################################################################################
#                                                     GSEA with multiple annotations
########################################################################################################################################

# geneList preparation

rm(all,gene_list,original_gene_list)

all <- as.data.frame(res_ouab_sh)
str(all)

# we want the log2 fold change 
original_gene_list <- all$log2FoldChange

# name the vector
names(original_gene_list) <- all$symbol

original_gene_list
# omit any NA values 
gene_list<-na.omit(original_gene_list)

anyDuplicated(names(gene_list))
gene_list <- gene_list[!duplicated(names(gene_list))]

gene_list

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

keytypes(org.Hs.eg.db)

################ GO BP

gseGOBP <- gseGO(geneList=gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont ="BP", 
                 keyType = "SYMBOL", 
                 nPerm = 100000, 
                 #minGSSize = 3, 
                 #maxGSSize = 1000,
                 verbose = TRUE,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.2,
                 by = "fgsea")

View(summary(gseGOBP))

dotplot(gseGOBP, showCategory=50, split=".sign") + facet_grid(.~.sign)

# GO:0010107 potassium ion import
# GO:1904064 positive regulation of cation transmembrane transport
# GO:0035794  positive regulation of mitochondrial membrane permeability
# GO:0008637  apoptotic mitochondrial changes
# GO:2001233  regulation of apoptotic signaling pathway
# GO:1901028  regulation of mitochondrial outer membrane permeabilization involved in apoptotic signalin
# GO:1902110	positive regulation of mitochondrial membrane permeability involved in apoptotic process
# GO:2001235	positive regulation of apoptotic signaling pathway
# GO:0006919	activation of cysteine-type endopeptidase activity involved in apoptotic process

gseaplot2(gseGOBP, 
          geneSetID = c("GO:0010107", "GO:1904064"), 
          color = c("#fa3c4c", "#d696bb"), 
          base_size = 14,
          #rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = TRUE,
          ES_geom = "line")

ggsave(file="./results_new/GSEA GO cations terms.jpeg", width=5.8, height=4.8, dpi=500)

# Just for visualization 
#gseGOBP2 <- gseGOBP
#View(gseGOBP2@result)
#rn <- rownames(gseGOBP2@result)
#rownames(gseGOBP2@result) <- c(1:891)
#gseGOBP2@result[216,2] <-"p1"
#gseGOBP2@result[73,2] <-"p2"
#gseGOBP2@result[237,2] <-"p5"
#gseGOBP2@result[385,2] <-"p4"
#gseGOBP2@result[540,2] <-"p3"
#rownames(gseGOBP2@result) <- rn

#gseaplot2(gseGOBP2,  #gseGOBP
#          geneSetID = c("GO:1902110", "GO:2001235", "GO:0006919"), 
#            color = c("#ffc100", "#ff7a7b", "#a24d53"),   
#          base_size = 14,
#          rel_heights = c(1.5, 0.5, 1),
#          #subplots = 1,
#          #pvalue_table = TRUE,
#          ES_geom = "line")

#ggsave(file="./results_new/GSEA GO apoptosis terms1.jpeg", width=5.8, height=4.8, dpi=500)

# for BI short presentation

gseaplot2(gseGOBP2,  #gseGOBP
          geneSetID = c("GO:1902110", "GO:0010107"), 
            color = c("#398564", "#f07f13"),   
          base_size = 14,
          rel_heights = c(1.5, 0.3, 1),
          #subplots = 1,
          #pvalue_table = TRUE,
          ES_geom = "line")

ggsave(file="./results_new/GSEA GO apoptosis terms6.jpeg", width=3.1, height=4.4, dpi=500)

# results export
gseGOBPs_df <- gseGOBP@result
str(gseGOBPs_df)
gseGOBPs_df<- arrange(gseGOBPs_df,p.adjust)
View(gseGOBPs_df)
write.csv(as.data.frame(gseGOBPs_df), file="./results_new/Supplementary.csv") # 

# additional visualization
dotplot(gseGOBP, showCategory=40, split=".sign") + facet_grid(.~.sign)
emapplot(gseGOBP)
cnetplot(gseGOBP, #categorySize="pvalue",
         showCategory = 2,
         foldChange=gene_list)
cnetplot(gseGOBP, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)
ridgeplot(gseGOBP)
heatplot(gseGOBP, foldChange=gene_list)
plotGOgraph(gseGOBP)
gseaplot(gseGOBPs, geneSetID = "GO:1904064", title = "positive regulation of cation transmembrane transport")

# potassium import related genes heatmap

# potassium import core genes - KCNJ2/WNK4/SLC12A8/WNK1/KCNJ8/KCNJ14/ABCC9/SLC12A7/SLC12A6
# pos reg cation transport - RELN/KCNJ2/WNK4/AGT/RGN/CX3CL1/PKD2/CAPN3/AMIGO1/TRPC1/HSPA2/WNK1/CACNA2D1

plotCounts(dds1, gene="ENSG00000123700.5", intgroup=c("CellType","Cells"), main = "KCNJ2 expression")
plotCounts(dds1, gene="ENSG00000163399.16", intgroup=c("CellType","Cells"), main = "ATP1A1 expression")
plotCounts(dds1, gene="ENSG00000143153.13", intgroup=c("CellType","Cells"), main = "ATP1B1 expression")
plotCounts(dds1, gene="ENSG00000129473.9", intgroup=c("CellType","Cells"), main = "ATP1B3 expression")

genes_hm1 <- c("ENSG00000123700.5","ENSG00000126562.17","ENSG00000221955.10", "ENSG00000060237.17",
        "ENSG00000121361.5", "ENSG00000182324.7", "ENSG00000069431.11", "ENSG00000113504.21", "ENSG00000140199.12",
        "ENSG00000189056.14", "ENSG00000135744.8", "ENSG00000130988.13",
        "ENSG00000006210.7", "ENSG00000118762.8", "ENSG00000092529.25", "ENSG00000181754.7", "ENSG00000144935.15",
        "ENSG00000126803.9", "ENSG00000153956.16")

rld_hm1 <- rld[genes_hm1,]
mat  <- assay(rld_hm1)
mat  <- mat - rowMeans(mat)
#anno <- as.data.frame(colData(rld_hm1)[, c("Cells","CellType","OuabSeno")])

rownames(anno) <- colnames(rld_hm1)
#colnames(anno)[1] <- "Cells"

colors<-colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(255) #length(breaksList)
anno <- as.data.frame(colData(rld_hm1)[, c("Cells","CellType","Senolysis")])
ann_colors <- list(Cells = c(Young = "#ffd596", Senescent = "#b2b2b2"),
                  CellType = c("END-MSCs" = "#b2d8d8", A549 = "#f7d0cb", "IMR-90" = "#e0ac69"),
                  Senolysis = c(Sensitive = "#ff556d", Resistant = "#9ace79"))

mat_eMSCs <- t(scale(t(mat[,1:8]))) 
mat_a549 <- t(scale(t(mat[,9:12])))
mat_imr <- t(scale(t(mat[,13:18])))
mat_scaled <- cbind(mat_eMSCs, mat_a549, mat_imr)
mat_scaled

plot5f <-pheatmap(mat_scaled,
         annotation_col = anno,
         #annotation_row = row_anno,
         #scale = "row",
         #breaks = breaksList,
         col=colors,
         annotation_colors = ann_colors,
         labels_row = rld_hm1@rowRanges$SYMBOL, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         gaps_col = c(8,12),
         show_colnames = FALSE,
         #legend_breaks = -1:4,
         main = "")

png(filename = "./results_new/HM1.png", width = 6200, height = 4400, units = "px", res = 1000)
plot5f
dev.off()

# apoptosis related genes heatmap

genes_hm2_df <- read.csv("apoptosis_hm.csv")
str(genes_hm2_df)
genes_hm2 <- levels(genes_hm2_df$symbol)
genes_hm2 <- rownames(subset(as.data.frame(res_ouab_sh), symbol %in% genes_hm2))

rld_hm2 <- rld[genes_hm2,]
mat  <- assay(rld_hm2)
mat  <- mat - rowMeans(mat)
#anno <- as.data.frame(colData(rld_hm1)[, c("Cells","CellType","OuabSeno")])

rownames(anno) <- colnames(rld_hm2)
#colnames(anno)[1] <- "Cells"

colors<-colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(255) #length(breaksList)
anno <- as.data.frame(colData(rld_hm1)[, c("Cells","CellType","Senolysis")])
ann_colors <- list(Cells = c(Young = "#ffd596", Senescent = "#b2b2b2"),
                   CellType = c("END-MSCs" = "#b2d8d8", A549 = "#f7d0cb", "IMR-90" = "#e0ac69"),
                   Senolysis = c(Sensitive = "#ff556d", Resistant = "#9ace79"))

mat_eMSCs <- t(scale(t(mat[,1:8]))) # apply(mat[,1:8], 2, scale) # scale(t(mat[,1:8])) apply(x, 1, function(x) x / sum(x, na.rm = TRUE))
mat_a549 <- t(scale(t(mat[,9:12])))
mat_imr <- t(scale(t(mat[,13:18])))
mat_scaled <- cbind(mat_eMSCs, mat_a549, mat_imr)
mat_scaled
# breaksList = seq(-2, 2, by = 1)

plot5g <-pheatmap(mat_scaled,
         annotation_col = anno,
         #annotation_row = row_anno,
         #scale = "row",
         #breaks = breaksList,
         col=colors,
         annotation_colors = ann_colors,
         labels_row = rld_hm2@rowRanges$SYMBOL, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         gaps_col = c(8,12),
         show_colnames = FALSE,
         #legend_breaks = -1:4,
         main = "")

# View(as.data.frame(rld_tf@rowRanges@elementMetadata))

png(filename = "./results_new/Figure 2f TFs HM.png", width = 6200, height = 14400, units = "px", res = 1000)
plot5g
dev.off()

############################################# Plotting KEGG Pathways
library(pathview)

# Plot specific KEGG pathways (with fold change) 
## pathway.id : KEGG pathway identifier
pathview(gene.data = gene_list, 
         pathway.id = "04218", 
         species = "hsa") # hsa04218 - Cellular senescence
pathview(gene.data = gene_matrix, 
         pathway.id = "04142", 
         species = "hsa") # hsa04142 - Lysosome
pathview(gene.data = gene_matrix, 
         pathway.id = "04020", 
         species = "hsa") # hsa04020 - Calcium signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04152", 
         species = "hsa") # ko04152 - AMPK signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04140", 
         species = "hsa") # hsa04140 - Autophagy
pathview(gene.data = gene_matrix, 
         pathway.id = "00010", 
         species = "hsa") # map00010 - Glycolysis / Gluconeogenesis
pathview(gene.data = gene_matrix, 
         pathway.id = "00020", 
         species = "hsa") # map00020 - Citrate cycle (TCA cycle)
pathview(gene.data = gene_matrix, 
         pathway.id = "03410", 
         species = "hsa") # ko03410 - Base excision repair
pathview(gene.data = gene_matrix, 
         pathway.id = "03450", 
         species = "hsa") # ko03450 - Non-homologous end-joining
pathview(gene.data = gene_matrix, 
         pathway.id = "02010", 
         species = "hsa") # ko02010 - ABC transporters
pathview(gene.data = gene_matrix, 
         pathway.id = "04330", 
         species = "hsa") # hsa04330 - Notch signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04310", 
         species = "hsa") # hsa04310 - Wnt signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04014", 
         species = "hsa") # hsa04014 - Ras signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04390", 
         species = "hsa") # hsa04390 - Hippo signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04350", 
         species = "hsa") # hsa04350 - TGF-beta signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04630", 
         species = "hsa") # hsa04630 - JAK-STAT signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04064", 
         species = "hsa") # hsa04064  - NF-kappa B signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04068", 
         species = "hsa") # hsa04068 - FoxO signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04151", 
         species = "hsa") # hsa04151 - PI3K-Akt signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04150", 
         species = "hsa") # hsa04150 - mTOR signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04512", 
         species = "hsa") # hsa04512 - ECM-receptor interaction
pathview(gene.data = gene_matrix, 
         pathway.id = "04514", 
         species = "hsa") # hsa04514 - Cell adhesion molecules
pathview(gene.data = gene_list, 
         pathway.id = "04210", 
         species = "hsa") # hsa04210 - Apoptosis
pathview(gene.data = gene_matrix, 
         pathway.id = "04110", 
         species = "hsa") # hsa04110 - Cell cycle
pathview(gene.data = gene_matrix, 
         pathway.id = "04115", 
         species = "hsa") # hsa04115 - p53 signaling pathway
pathview(gene.data = gene_matrix, 
         pathway.id = "04550", 
         species = "hsa") # hsa04550 - Signaling pathways regulating pluripotency of stem cells
pathview(gene.data = gene_list, 
         pathway.id = "04978", 
         species = "hsa") # hsa04978 - Mineral absorption
pathview(gene.data = gene_matrix, 
         pathway.id = "04750", 
         species = "hsa") # hsa04750 - Inflammatory mediator regulation of TRP channels
pathview(gene.data = gene_matrix, 
         pathway.id = "04136", 
         species = "hsa") # hsa04136 - Autophagy - other


