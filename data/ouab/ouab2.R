
###########################################################################################################################################

#                                                   RNA-seq analysis for OUABAIN THEME

###########################################################################################################################################

############## Downloading and preprocessing of raw data from the SRA 

prefetch SRR5931953  # A549 - ctr - GSE102639
         SRR5931954  # A549 - ctr - GSE102639 
         SRR5931961  # A549 - sen etoposide - GSE102639 
         SRR5931962  # A549 - sen etoposide - GSE102639 

         SRR8145428 # IMR-90 - ctr - GSE122081
         SRR8145429 # IMR-90 - ctr - GSE122081 
         SRR8145430 # IMR-90 - ctr - GSE122081 
         SRR8145446 # IMR-90 - sen OIS - GSE122081
         SRR8145447 # IMR-90 - sen OIS - GSE122081 
         SRR8145448 # IMR-90 - sen OIS - GSE122081

find . -maxdepth 2 -name "*.sra" -exec cp {} ~/mks/ouab_theme/raw_reads \;

fasterq-dump --split-files SRR5931953.sra # convertation to .fastq
rm *.sra # removing .sra files
gzip # each .fastq file

fastqc *.gz # raw data QC
find . -maxdepth 2 -name "*.html" -exec mv {} ~/mks/ouab_theme/step_1_fastqc_raw_reads \; 
find . -maxdepth 2 -name "*.zip" -exec mv {} ~/mks/ouab_theme/step_1_fastqc_raw_reads \;

# Additional trimming of "bad" IMR-90 dataset using bbduk

# ~/Software/BBMap_38.75/bbmap/bbduk.sh in1=SRR9016157.sra_1.fastq.gz in2=SRR9016157.sra_2.fastq.gz out1=SRR9016157.1.fastq.gz out2=SRR9016157.2.fastq.gz ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ftl=10

~/Software/BBMap_38.75/bbmap/bbduk.sh in=SRR8145428.fastq.gz out=SRR8145428.fastq.gz ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ftl=10

################## transcriptome-based "aligmnent" with salmon

# for pair-end reads

salmon quant -i ~/mks/data_and_results/step9_salmon_quant/salmon_index_sa_our/salmon_index_sa \
-l A \
-1 <(zcat SRR9016157.1.fastq.gz) \
-2 <(zcat SRR9016157.2.fastq.gz) \
-p 12 \
-o ~/mks/ouab_theme/step_2_salmon/quants/SRR9016157 \
--numBootstraps 30 \
--seqBias \
--gcBias \
--validateMappings

# or by the script, for single-end data: (remember to chmod +x)

#!/bin/bash

for data in *;
do
if [[ ${data} == *.fastq.gz ]]
then
echo "Processing sample ${data}"
salmon quant -i ~/mks/data_and_results/step9_salmon_quant/salmon_index_sa_our/salmon_index_sa \
-l A \
-r <(zcat ${data}) \
-p 12 \
-o ~/mks/ouab_theme/step_2_salmon/quants/$(echo ${data} | cut -c 1-10) \
--numBootstraps 30 \
--seqBias \
--gcBias \
--validateMappings
else
  echo "File is not reads data"
fi
done

#######################  Tximeta quantification 

setwd("~/mks/ouab_theme/step_2_salmon") # setting wd with the folder "quants" and predesigned "coldata.csv" file
library(Rmisc)
library(xlsx)
library(ggplot2)
library(dplyr)

dir <- ("quants")
list.files(file.path(dir))

dir2 <- (".")
csvfile <- file.path(dir2, "coldatadf3.csv")
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)

str(coldata)
coldata$SampleID <- as.character(coldata$SampleID)
coldata$names <- coldata$SampleGName 
coldata$files <- file.path(dir2, "quants", coldata$names, "quant.sf")
file.exists(coldata$files)
coldata # look at the data

library(tximeta)

se <- tximeta(coldata) # quantification to transcript-base level
dim(se) # look at the munber of dimentions
head(rownames(se))

gse <- summarizeToGene(se) # summarization to gene-level

library(ggplot2)
library(dplyr)

library(org.Hs.eg.db)
library(AnnotationDbi)

keytypes(org.Hs.eg.db) # add nessessary annotations
gse <- addIds(gse, "REFSEQ", gene=TRUE)
gse <- addIds(gse, "ENTREZID", gene=TRUE)
gse <- addIds(gse, "GENENAME", gene=TRUE)
gse <- addIds(gse, "SYMBOL", gene=TRUE)

mcols(gse)

######################## EVA and DE analysis with DESeq2

library(DESeq2)

### Construction of the DESeqDataSet - the DESeq2 custom class object (almost* the same as Bioconductor custom class object SummarizedExperiment)

str(gse@colData$CellType) 
gse@colData$CellType <- as.factor(gse@colData$CellType)
gse@colData$Cells <- as.factor(gse@colData$Cells)
gse@colData$OuabSeno <- as.factor(gse@colData$OuabSeno)
levels(gse@colData$CellType)
levels(gse@colData$Cells)
levels(gse@colData$OuabSeno)
colData(gse)$CellType <- relevel(colData(gse)$CellType, ref = "eMSCs")
levels(gse@colData$CellType)

dds <- DESeqDataSet(gse, design = ~ Cells + OuabSeno + Cells:OuabSeno)
dds # get summary
nrow(dds) # learn total number of rows

### EVA 

dds1 <- dds[ rowSums(counts(dds)) > 1, ] # removing rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds1)
# or
keep <- rowSums(counts(dds) >= 5) >= 4 # minimal filtering to reduce the size of the dataset. We do not need to retain genes 
# if they do not have a count of 5 or more for 4 or more samples 
# as these genes will have no statistical power to detect differences, 
# and no information to compute distances between samples
table(keep)
dds1 <- dds[keep,]

rld <- DESeq2::rlog(dds1,blind = FALSE) # DESeq2 transformation for count data that stabilize the variance across the mean - the regularized logarithm (rlog)
# we specified blind = FALSE, which means that differences between cell lines and senescence (the variables in the design) 
# will not contribute to the expected variance-mean trend of the experiment. 
# The experimental design is not used in the transformation directly* 
class(rld) # DESeqTransform object which is based on the SummarizedExperiment class (NOT SUITABLE FOR DE ANALYSIS!)

DESeq2::plotPCA(rld, intgroup = c("CellType", "Cells"))

# Senescence validation
total <- as.data.frame(dds1@rowRanges)
View(total)
total <- total[,c(6,8)]
total1 <- total[,c(6,10)]
sub1 <- read.xlsx("Cellular senescence keggmap=04218.xlsx", sheetIndex = 1)
str(sub1)
str(total)
sub1$ENTREZID <- as.character(sub1$ENTREZID)
merged <- merge(total,sub1,by="ENTREZID")
head(merged)
#merged$ensembl <- as.character(merged$ensembl)
str(merged)
sub1v <- merged[,2]

rld_sen <- rld[sub1v,]

DESeq2::plotPCA(rld_sen, intgroup = c("CellType", "Cells"))

############ DEG analysis with DESeq2

dds <- DESeq(dds) # running DE analysis (on non-transformed/normilized counts!!!)

resultsNames(dds) # get terms from the model

res_ouab <- results(dds, name = "Cellssen.OuabSenonot_affected") # DEGs for interaction

# for now the difference reflects lines vs eMSCs -> *(-1) lfc

mcols(res_ouab, use.names = TRUE) # description of DataFrames
summary(res_ouab) # main statistics

plotMA(res_ouab, ylim = c(-10, 10)) # look at all genes
#abline(h=c(-1,1), col="dodgerblue", lwd=2)

library("ashr") # advanced lcf shrinkage method (more informative visualization and more accurate ranking of genes by effect size - adequate lcf values)
res_ouab_sh <- lfcShrink(dds, coef = 4, type="apeglm") 

plotMA(res_ouab_sh, ylim = c(-10, 10))
#abline(h=c(-1,1), col="dodgerblue", lwd=2)

### Annotations with different IDs

library("AnnotationDbi")
library("org.Hs.eg.db")
library("GenomicFeatures")

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

res_ouab_shs <- subset(res_ouab_sh, padj < 0.1) # subset of the BH significant (FDR=0.1)
summary(res_ouab_shs)

res_ouab_shs <- res_ouab_shs[order(res_ouab_shs$log2FoldChange),]
View(as.data.frame(res_ouab_shs)) # look manually for interesting genes
plotCounts(dds, gene="ENSG00000157445.15", intgroup=c("CellType","Cells"))

############### subsetting 

setwd("~/mks/ouab_theme/step_2_salmon/results for ouab")

res_ouab_shs <- res_ouab_shs[order(res_ouab_shs$log2FoldChange),] # order by lfc
write.csv(as.data.frame(res_ouab_shs), file="All BH-sig DEGs for ouab theme.csv") # ALL BH-sig DEGs, ordered by lfc

res_ouab_shs_up <- subset(res_ouab_shs, log2FoldChange > 0) # subset of the BH significant (FDR=0.1) genes up-regulated at least by 1,5 times
res_ouab_shs_up <- res_ouab_shs_up[order(res_ouab_shs_up$log2FoldChange, decreasing = TRUE), ] # order by lfc
res_ouab_shs_dw <- subset(res_ouab_shs, log2FoldChange < 0) # subset of the BH significant (FDR=0.1) genes down-regulated at least by 1,5 times
res_ouab_shs_dw <- res_ouab_shs_dw[order(res_ouab_shs_dw$log2FoldChange),] # order by lfc

write.csv(as.data.frame(res_ouab_shs_up), file="Up-reg BH-sig DEGs for ouab theme.csv") 
write.csv(as.data.frame(res_ouab_shs_dw), file="Down-reg BH-sig DEGs for ouab theme.csv") 

########################################################################################################################################
#                                                     GSEA with multiple annotations
########################################################################################################################################

library(Rmisc)
library(xlsx)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

#################### geneList preparation

rm(all,gene_list,original_gene_list)

all <- as.data.frame(res_ouab_shs)
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
library(DOSE)
library(topGO)

################ GO BP

gseGOBP <- gseGO(geneList=gene_list,
             OrgDb = org.Hs.eg.db,
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 1000,
             verbose = TRUE,
             pAdjustMethod = "BH",
             pvalueCutoff = 1,
             by = "fgsea")

View(summary(gseGOBP))

gseGOBPs <- simplify(gseGOBP,
                  cutoff = 0.7,
                  by = "p.adjust",
                  select_fun = min,
                  measure = "Wang",
                  semData = NULL)

#GO:1904064 

View(summary(gseGOBPs2))

dotplot(gseGOBP, showCategory=50, split=".sign") + facet_grid(.~.sign)
dotplot(gseGOBPs2, showCategory=30, split=".sign") + facet_grid(.~.sign)

# results export
gseGOBPs_df <- gseGOBP@result
str(gseGOBPs_df)
gseGOBPs_df<- arrange(gseGOBPs_df,p.adjust)
View(gseGOBPs_df)
write.xlsx(gseGOBPs_df,"gseGOBP.xlsx")

# additional visualization
dotplot(gse, showCategory=40, split=".sign") + facet_grid(.~.sign)
emapplot(gseGOBP)
cnetplot(gseGOBP, #categorySize="pvalue", 
         foldChange=gene_list)
cnetplot(gseGOBP, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)
ridgeplot(gseGOBP)
heatplot(gseGOBP, foldChange=gene_list)
plotGOgraph(gseGOBP)
    gseaplot(gseGOBPs, geneSetID = "GO:1904064", title = "positive regulation of cation transmembrane transport")
#gseaplot2(gse, geneSetID = "GO:0043065", title = gse$Description [ID = "GO:0043065"])

# GO CC
gseGOCC <- gseGO(geneList=gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont ="CC", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 1000,
                 verbose = TRUE,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 by = "fgsea")

View(summary(gseGOCC))

gseGOCCs <- simplify(gseGOCC,
                     cutoff = 0.7,
                     by = "p.adjust",
                     select_fun = min,
                     measure = "Wang",
                     semData = NULL)

View(summary(gseGOCCs))

dotplot(gseGOCC, showCategory=50, split=".sign") + facet_grid(.~.sign)
dotplot(gseGOBPs, showCategory=30, split=".sign") + facet_grid(.~.sign)

# results export
gseGOCC_df <- gseGOCC@result
gseGOCC_df<- arrange(gseGOCC_df,p.adjust)
View(gseGOCC_df)
write.xlsx(gseGOCC_df,"gseGOCC.xlsx")

# additional visualization
dotplot(gse, showCategory=40, split=".sign") + facet_grid(.~.sign)
emapplot(gseGOBPs)
cnetplot(gse, #categorySize="pvalue", 
         foldChange=gene_list)
cnetplot(gse, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)
ridgeplot(gse)
heatplot(gse, foldChange=gene_list)
plotGOgraph(gse)
gseaplot(gse, geneSetID = "GO:0043065")
gseaplot2(gse, geneSetID = "GO:0043065", title = gse$Description [ID = "GO:0043065"])


# GO MF
gseGOMF <- gseGO(geneList=gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont ="MF", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 1000,
                 verbose = TRUE,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 by = "fgsea")

View(summary(gseGOMF))

gseGOMFs <- simplify(gseGOMF,
                     cutoff = 0.7,
                     by = "p.adjust",
                     select_fun = min,
                     measure = "Wang",
                     semData = NULL)

View(summary(gseGOMFs))

dotplot(gseGOMF, showCategory=50, split=".sign") + facet_grid(.~.sign)
dotplot(gseGOBPs, showCategory=30, split=".sign") + facet_grid(.~.sign)

# results export
gseGOMF_df <- gseGOMF@result
gseGOMF_df<- arrange(gseGOMF_df,pvalue)
View(gseGOMF_df)
write.xlsx(gseGOMF_df,"gseGOMF.xlsx")

# additional visualization
dotplot(gse, showCategory=40, split=".sign") + facet_grid(.~.sign)
emapplot(gseGOMF)
cnetplot(gse, #categorySize="pvalue", 
         foldChange=gene_list)
cnetplot(gse, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)
ridgeplot(gse)
heatplot(gse, foldChange=gene_list)
plotGOgraph(gse)
gseaplot(gse, geneSetID = "GO:0043065")
gseaplot2(gse, geneSetID = "GO:0043065", title = gse$Description [ID = "GO:0043065"])

# KEGG 

rm(all,gene_list,original_gene_list)

all <- as.data.frame(res_ouab_shs)
str(all)

# we want the log2 fold change 
original_gene_list <- all$log2FoldChange

# name the vector
names(original_gene_list) <- all$entrez

original_gene_list
# omit any NA values 
gene_list<-na.omit(original_gene_list)

anyDuplicated(names(gene_list))
gene_list <- gene_list[!duplicated(names(gene_list))]

gene_list

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = gene_list,
               organism     = 'hsa',
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 1000,
               verbose = TRUE,
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               by = "fgsea")

View(summary(kk2))
dotplot(kk2, showCategory=50, split=".sign") + facet_grid(.~.sign)
emapplot(kk2)

kk2r <- setReadable(kk2, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
View(summary(kk2r))

# results export
kk2r_df <- kk2r@result
kk2r_df<- arrange(kk2r_df,pvalue)
View(kk2r_df)
write.xlsx(kk2r_df,"KEGG.xlsx")

# KEGG Module Gene Set Enrichment Analysis

mkk2 <- gseMKEGG(geneList = gene_list,
                 organism = 'hsa',
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 1000,
                 verbose = TRUE,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 by = "fgsea")

View(summary(mkk2))
dotplot(mkk2, showCategory=50, split=".sign") + facet_grid(.~.sign)
emapplot(mkk2)

mkk2r <- setReadable(mkk2, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
View(summary(mkk2r))

# results export
mkk2r_df <- mkk2r@result
mkk2r_df<- arrange(mkk2r_df,pvalue)
View(mkk2r_df)
write.xlsx(mkk2r_df,"MKEGG.xlsx")

# Reactome

library(ReactomePA)
react <- gsePathway(geneList = gene_list,
                     organism = 'human',
                     nPerm = 10000, 
                     minGSSize = 3, 
                     maxGSSize = 1000,
                     verbose = TRUE,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     by = "fgsea")
View(summary(react))

reactr <- setReadable(react, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
View(summary(reactr))

dotplot(reactr, showCategory=50, split=".sign") + facet_grid(.~.sign)

# results export
reactr_df <- reactr@result
reactr_df<- arrange(reactr_df,pvalue)
View(reactr_df)
write.xlsx(reactr_df,"Reactome.xlsx")

# WP

wpgmtfile <- system.file("extdata/wikipathways-20200810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

gseWP <- GSEA(gene_list, 
             TERM2GENE = wpid2gene, 
             TERM2NAME = wpid2name, 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 1000,
             verbose = TRUE,
             pAdjustMethod = "BH",
             pvalueCutoff = 1,
             by = "fgsea")
head(gseWP)

View(summary(gseWP))

gseWP <- setReadable(gseWP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
View(summary(gseWP))

dotplot(gseWP, showCategory=50, split=".sign") + facet_grid(.~.sign)

# results export
gseWP_df <- gseWP@result
gseWP_df<- arrange(gseWP_df,pvalue)
View(gseWP_df)
write.xlsx(gseWP_df,"WikiPathways.xlsx")


# MSigDB
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msigdbr")

library(msigdbr)
msigdbr_show_species()

m_df <- msigdbr(species = "Homo sapiens")%>% 
  dplyr::select(gs_name, entrez_gene)
head(m_df)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

msigdb <- GSEA(geneList = gene_list, 
            TERM2GENE = m_df,
            nPerm = 10000, 
            minGSSize = 3, 
            maxGSSize = 1000,
            verbose = TRUE,
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            by = "fgsea")

View(summary(msigdb))

msigdbr <- setReadable(msigdb, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
View(summary(msigdbc6r))

dotplot(msigdbr, showCategory=50, split=".sign") + facet_grid(.~.sign)

# results export
msigdbr_df <- msigdbr@result
msigdbr_df<- arrange(msigdbr_df,pvalue)
View(msigdbr_df)
write.xlsx(msigdbr_df,"MSigBP.xlsx")


########################################################################################################################################
#                                                 Specific visualization based on transformed data
########################################################################################################################################

rld <- DESeq2::rlog(dds,blind = FALSE) # all data

plotCounts(dds, gene=which.min(res_ouab_s$padj), intgroup=c("CellType","Cells")) # validation of dds design

library(pheatmap)

mat = assay(rld)[ head(order(res_ouab_s$padj),5), ] # select the top 5 genes with the lowest padj
anno <- as.data.frame(colData(rld)[, c("Cells","CellType","Ouab")])
pheatmap(mat, 
         annotation_col=anno,
         scale = "row",
         labels_row = rld@rowRanges$GENENAME, 
         cluster_cols = F,
         gaps_col = c(8,12),
         show_colnames = FALSE,
         main = "")

plotCounts(dds, gene="ENSG00000280339.1", intgroup=c("CellType","Cells"))
#= dds$Description [ID = "GO:0043065"]

################################### Making subset from DEGs for GO relevant genes

total <- as.data.frame(rld@rowRanges)
total2 <- total[,c(6,8,10)]
View(total2)

sub3 <- read.xlsx("uniprot_GO:1904064.xlsx", sheetIndex = 1)
str(sub3)
str(total2)
sub3$ENTREZID <- as.character(sub3$ENTREZID)
merged <- merge(total2,sub3,by="ENTREZID")
head(merged)
merged_ions <- merged$gene_id

rld_ions <- rld[merged_ions,]

mat  <- assay(rld_ions)
anno <- as.data.frame(colData(rld_ions)[, c("Cells","CellType","OuabSeno")])

pheatmap(mat,
         annotation_col = anno,
         scale = "row",
         labels_row = rld_ions@rowRanges$GENENAME, 
         cluster_cols = F,
         gaps_col = c(8,12),
         show_colnames = FALSE, 
         main = "")




############################################# Plotting KEGG Pathways
library(pathview)

# Plot specific KEGG pathways (with fold change) 
## pathway.id : KEGG pathway identifier
pathview(gene.data = gene_matrix, 
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
pathview(gene.data = gene_matrix, 
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
pathview(gene.data = gene_matrix, 
         pathway.id = "04978", 
         species = "hsa") # hsa04978 - Mineral absorption
pathview(gene.data = gene_matrix, 
         pathway.id = "04750", 
         species = "hsa") # hsa04750 - Inflammatory mediator regulation of TRP channels
pathview(gene.data = gene_matrix, 
         pathway.id = "04136", 
         species = "hsa") # hsa04136 - Autophagy - other


