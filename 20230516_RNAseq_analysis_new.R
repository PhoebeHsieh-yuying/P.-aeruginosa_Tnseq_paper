##NOTE: 20230516 when I made a heatmap of candidate genes from Tn-seq, I found my dds file missed a couple genes. Later I found the gtf and gff files downloaded from Pseudomonas.com contain different numbers of annotated genes. To solve this, Janet helped me to use "20230516_checkGFF_GTF_files.R" to create a correct GTF file for running featurecounts (sbatch RNAseq_featurecounts_new.sh). I will get new output files from featurecounts with locus(PAXXXX) as row names, then I can run the below R lines to read all 12 samples and combine them into one big datafile.


## Load Libraries
library(tidyverse)
library(DESeq2)

## clear environment
rm(list=ls())

## define some files and locations
# my mac does not mount the shared drive in a consistent location - annoying! 
# so I have first figure out where it is:
sharedDriveLocation <- "/Volumes/fh/fast/malik_h/"
if (!file.exists(sharedDriveLocation)) { sharedDriveLocation <- "/Volumes/malik_h/" }

## change working dir
workingDir <- paste(sharedDriveLocation, "/user/yhsieh/Raw_Sequence_Reads/RNA_seq_genewiz_2023_0201/DESeq2_analysis_new", sep="")
setwd(workingDir)

## specify where Phoebe's analysis files are 
inputFileDir <- paste(sharedDriveLocation, "", sep="/user/yhsieh/Raw_Sequence_Reads/RNA_seq_genewiz_2023_0201/DESeq2_analysis_new")

## read the fixed output of featureCounts. Use lapply to read each file, but lapply produces a list, which is not good for downstream analysis(for that, use sapply instead)
dataFiles <- list.files(inputFileDir, pattern ="_featureCounts_newoutput.txt", full.names = TRUE)

dataList <- lapply(dataFiles, function(x) {
  Y <- read.delim(file = x,
                  header = TRUE, # If you have a header row you may want to turn this to true, and remove the "skip" param
                  skip = 1, # there are additional headers in these files we want to skip - in ours there are 4 rows to skip but if your header is the first row, you can remove this param or set it to 0
                  sep = "\t",
                  stringsAsFactors = FALSE)
  Y
})


## extract sample name from each file in dataList. information of sample name is in column 7, and then use gsub to remove unwanted strings (_subread_sorted.bam)
sampleNames <- sapply(dataList, function(x) { colnames(x)[7] }) 

sampleNames <- gsub("_subread_sorted.bam", "", sampleNames)

## counts_matrix is a combined table with counts from 12 samples and information of Geneid (from the first file in datalist). new output files have PAXXXX as Geneid
counts_matrix <- sapply( dataList, function(x) {
  x[,7]
} )
colnames(counts_matrix) <- sampleNames
rownames(counts_matrix) <- dataList[[1]][,"Geneid"]

# 20240515 write csv file of counts_matrix for raw data of Table S1
write.csv (counts_matrix, file = "raw_data_TableS1.csv")

### the below script can allow me to get function for each genes of interest, based on locus_tag. I can use this script to annotate the function of a list of genes after DESeq2. I need to change the below script as my new output files don't have annotation for gene's function
fakeResults <- as.data.frame(counts_matrix[4000:4006,1:3])

fakeResults[,"locus_tag"] <- rownames(fakeResults)

left_join( fakeResults, dataList[[1]][,c("locus_tag","function.")], by="locus_tag")


## Make a dataframe for design table of DESeq2(similar as coldata). After getting this design table (coldata), I can run DESeq2  using counts_matrix and coldata

## test if a row has Y1, if yes, then TRUE; overwise, FALSE. as.character converts logic values into characters. culture_conditions and Mg_conditions are the contents for two columns in coldata
culture_conditions <- as.character(grepl("Y1", sampleNames)) 
culture_conditions <- gsub("FALSE", "monoculture",culture_conditions)
culture_conditions <- gsub("TRUE", "coculture",culture_conditions)
#culture_conditions <- factor(culture_conditions)
culture_conditions <- factor(culture_conditions, levels=c("monoculture", "coculture"))

Mg_conditions <- as.character(grepl("Mg", sampleNames))
Mg_conditions <- gsub("FALSE", "no_Mg",Mg_conditions)
Mg_conditions <- gsub("TRUE", "with_Mg",Mg_conditions)
Mg_conditions <- factor(Mg_conditions, levels =c("no_Mg", "with_Mg"))

treatments <- factor(c("BHI", "BHI", "BHI", "BHIMg", "BHIMg", "BHIMg", "BHIY1", "BHIY1", "BHIY1", "BHIY1Mg", "BHIY1Mg", "BHIY1Mg"), levels = c("BHI", "BHIMg", "BHIY1", "BHIY1Mg"))

coldata <- data.frame(sampleName=sampleNames,
                      #get last character of sampleNames. create the content of replicate
                      replicate= sapply(sampleNames, function(x) {
                        len_x <- nchar(x)
                        rep <- substr(x, start=len_x, stop=len_x)
                        rep <- as.integer(rep)
                        return(rep)
                      }),
                      culture_condition = culture_conditions,
                      Mg_condition = Mg_conditions,
                      treatment = treatments
                      
)


write.csv(coldata, file ="designtable_RNAseq_20230505.csv")

### Run DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~culture_condition + Mg_condition)

##20230510, the belwo script is suggested by Janet. I need to specify dds for each trial
# the second fator in design is the primary parameter to compare by default
# dds_twoThingsNoInteraction <- DESeqDataSetFromMatrix(countData = counts_matrix,
#                                                      colData = coldata,
#                                                      design = ~culture_condition + Mg_condition) 
# dds_twoThingsNoInteraction <- DESeq(dds_twoThingsNoInteraction)
# 
# dds_groups <- DESeqDataSetFromMatrix(countData = counts_matrix,
#                                      colData = coldata,
#                                      design = ~ group) 
# dds_groups <- DESeq(dds_groups)
# 
# dds_twoThingsWithInteraction <- DESeqDataSetFromMatrix(
#   countData = counts_matrix,
#   colData = coldata,
#   design = ~culture_condition + Mg_condition + culture_condition:Mg_condition)
# dds_twoThingsWithInteraction <- DESeq(dds_twoThingsWithInteraction)
# res_twoThingsWithInteraction_cult <- results(dds_twoThingsWithInteraction,
#                                              name="culture_condition_coculture_vs_monoculture")
# 
# res_twoThingsWithInteraction_Mg <- results(dds_twoThingsWithInteraction,
#                                            name="Mg_condition_with_Mg_vs_no_Mg")
# 
# res_twoThingsWithInteraction_interaction <- results(dds_twoThingsWithInteraction,
#                                                     name="culture_conditioncoculture.Mg_conditionwith_Mg")
# 
# 
# interactivate(dds_twoThingsWithInteraction)

#res_twoThingsWithInteraction_cult["PA4825",]
#res_twoThingsWithInteraction_Mg["PA4825",]
#res_twoThingsWithInteraction_interaction["PA4825",]
#dds <- DESeq(dds)
#resultsNames(dds)

## interaction terms to the design formula to compare DEG between co-culture and monoculture in BHI (or BHI+Mg). This is to check the effect of co-culture for different sets of samples(BHI, BHIMg). I follow the paragraph of "interactions" in the DESeq2 document
dds$group <- factor(paste0(dds$culture_condition, dds$Mg_condition))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
res <-results(dds)
## Pre-filtering. Remove rows that have reads smaller than 10 per sample. I have 12 samples, 
keep <- rowSums(counts(dds)) >= 120
dds <- dds[keep,]

## Find DEGs between BHIY1 and BHI
resBHI <- results(dds, contrast=c("group", "cocultureno_Mg", "monocultureno_Mg"))
resBHI_FLT <- resBHI[which(abs(resBHI[, "log2FoldChange"]) >2& resBHI[,"padj"] < 0.05),]
resBHI_FLT
summary(resBHI_FLT)
resBHI_FLT <- as.data.frame(resBHI_FLT)

#Make a volcano plot between BHIY1 and BHI, with labeling PA4824, PA4825, and PA4826
library(EnhancedVolcano)

keyvals2 <- ifelse(
  resBHI$log2FoldChange < -2 &resBHI$padj<0.05, 'royalblue',
  ifelse(resBHI$log2FoldChange > 2 &resBHI$padj
         <0.05, '#A60021',
         '#E5E5E5'))
         keyvals2[is.na(keyvals2)] <- '#E5E5E5'
           names(keyvals2)[keyvals2 == '#A60021'] <- 'Upregulated genes'
           names(keyvals2)[keyvals2 == '#E5E5E5'] <- 'not significant'
           names(keyvals2)[keyvals2 == 'royalblue'] <- 'Downregulated genes'
           
           p1 <- EnhancedVolcano(resBHI,
                                 lab = rownames(resBHI),
                                 x = 'log2FoldChange',
                                 y = 'pvalue',
                                 selectLab =c('PA4825','PA4824','PA4826'),
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 #ylim =c(0,180), xlim = c(-10,25),
                                 title = 'Gene Expression in Co-culture versus Monoculture in BHI',
                                 pCutoff = 0.05, pCutoffCol = 'padj',
                                 FCcutoff = 2.0,
                                 pointSize = 4,
                                 labSize = 4,
                                 colCustom = keyvals2,
                                 colAlpha = 1,
                                 subtitle = paste0('p-value cutoff (black dotted line) drawn ',
                                                   'at equivalent of adjusted p=0.05'),
                                 #legendLabSize = 15,
                                 #legendPosition = 'left',
                                 legendIconSize = 5.0,
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.5,
                                 colConnectors = 'grey50',
                                 gridlines.major = TRUE,
                                 gridlines.minor = FALSE,
                                 border = 'partial',
                                 borderWidth = 1.5,
                                 borderColour = 'black')
           p1

# create an annotation file for matching Geneid with gune function
gff_annotation <- read.csv ("PAO1_gff_annotation.csv", header = TRUE)
resBHI_FLT[,"Geneid"] <- rownames(resBHI_FLT)
resBHI_FLT_join <- left_join( resBHI_FLT, gff_annotation[,c("Geneid","function.")], by="Geneid")

write.csv(resBHI_FLT_join,file= "RNAseq_DEgenes_BHIY1_BHI_4FC_padj0.05_filter120_function.csv") 

# create an annotation file of each gene with FC and p-value. I need to convert resBHI into a dataframe before getting a list of genes
gff_annotation <- read.csv ("PAO1_gff_annotation.csv", header = TRUE)
resBHI[,"Geneid"] <- rownames(resBHI)
resBHI <- as.data.frame(resBHI)

resBHI_join <- left_join( resBHI, gff_annotation[,c("Geneid","function.")], by="Geneid")

write.csv(resBHI_join,file= "RNAseq_DEgenes_BHIY1_BHI_4FC_filter120_function.csv") 

## Find DEGs between MgY1 and Mg
resMg <- results(dds, contrast=c("group", "coculturewith_Mg", "monoculturewith_Mg"))
resMg_FLT <- resMg[which(abs(resMg[, "log2FoldChange"]) >2& resMg[,"padj"] < 0.05),]
resMg_FLT
summary(resMg_FLT)

resMg_FLT <- as.data.frame(resMg_FLT)

resMg_FLT[,"Geneid"] <- rownames(resMg_FLT)
resMg_FLT_join <- left_join( resMg_FLT, gff_annotation[,c("Geneid","function.")], by="Geneid")

write.csv(resMg_FLT_join,file= "RNAseq_DEgenes_MgY1_Mg_4FC_padj0.05_filter120_function.csv") 

sum(resBHI$padj < 0.05, na.rm =TRUE)
summary(resMg)


## Find DEGs between co-culture and monoculture
dds$culture_condition <- relevel(dds$culture_condition, ref = "monoculture")
dds <- DESeq(dds)
resCO <- results(dds, contrast= c("culture_condition", "coculture", "monoculture"))
resCO_FLT <- resCO[which(abs(resCO[, "log2FoldChange"]) >2& resCO[,"padj"] < 0.05),]
resCO_FLT

resCO_FLT <- as.data.frame(resCO_FLT)

resCO_FLT[,"locus_tag"] <- rownames(resCO_FLT)
resCO_FLT_join <- left_join( resCO_FLT, dataList[[1]][,c("locus_tag","function.")], by="locus_tag")

write.csv(resCO_FLT_join,file= "RNAseq_DEgenes_CO_MONO_4FC_padj0.05_filter120_function.csv")

## plot MA plot 
plotMA(resBHI)

plotMA(resMg)

# make PCA plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("culture_condition", "Mg_condition"))

## plot heatmap of dds(counts_matirx)
ntd <- normTransform(dds)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("culture_condition","Mg_condition")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$culture_condition, vsd$Mg_condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


## Make heatmap of genes of interest across conditition. I should use the following script for making heatmaps for publication
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("genefilter")
BiocManager::install("gplots")
BiocManager::install("RcolorBrewer")

# for dds using group design
rld <- rlog(dds, blind=FALSE)
topVarGenes<-head(order(rowVars(assay(rld)),decreasing=TRUE),60)

rld_twoThingsWithInteraction <- rlog(dds_twoThingsWithInteraction, blind=FALSE)
topVarGenes_twoThingsWithInteraction <- head(order(rowVars(assay(rld_twoThingsWithInteraction)),decreasing=TRUE),200)


library("genefilter") 
library("gplots")
library("RColorBrewer")


heatmap.2(assay(rld)[topVarGenes,],
          scale="row", trace="none",
          dendrogram="both" , 
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

heatmap.2(assay(rld)[topVarGenes,],
          scale="row", trace="none",
          Rowv = TRUE,
          Colv =TRUE,
          distfun = dist,
          hclustfun = hclust,
          dendrogram="row" , 
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

### top 200 most variable genes
heatmap.2(assay(rld_twoThingsWithInteraction)[topVarGenes_twoThingsWithInteraction,],
          scale="row", trace="none",
          dendrogram="both" , 
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

### get gene lists

# table( abs(res_twoThingsWithInteraction_cult$log2FoldChange)>=1 )
# table( res_twoThingsWithInteraction_cult$padj < 0.05 )
# 
# table( abs(res_twoThingsWithInteraction_cult$log2FoldChange)>=2 )
# 
# table( abs(res_twoThingsWithInteraction_cult$log2FoldChange)>=2 & res_twoThingsWithInteraction_cult$padj < 0.05 )
# 
# 
# res_twoThingsWithInteraction_cult_interestingGenes <- res_twoThingsWithInteraction_cult[which(abs(res_twoThingsWithInteraction_cult[, "log2FoldChange"]) >2 & res_twoThingsWithInteraction_cult[,"padj"] < 0.05),]
# 
# res_twoThingsWithInteraction_Mg_interestingGenes <- res_twoThingsWithInteraction_Mg[which(abs(res_twoThingsWithInteraction_Mg[, "log2FoldChange"]) >1 & res_twoThingsWithInteraction_Mg[,"padj"] < 0.05),]
# 
# res_twoThingsWithInteraction_interaction_interestingGenes <- res_twoThingsWithInteraction_interaction[which(abs(res_twoThingsWithInteraction_interaction[, "log2FoldChange"]) >1 & res_twoThingsWithInteraction_interaction[,"padj"] < 0.05),]
# 
# res_twoThingsWithInteraction_interestingGenes_allComparison <- c(rownames(res_twoThingsWithInteraction_cult_interestingGenes), rownames(res_twoThingsWithInteraction_Mg_interestingGenes), rownames(res_twoThingsWithInteraction_interaction_interestingGenes))
# 
# res_twoThingsWithInteraction_interestingGenes_allComparison <- unique(res_twoThingsWithInteraction_interestingGenes_allComparison)


resBHI <- results(dds, contrast=c("group", "cocultureno_Mg", "monocultureno_Mg"))
resBHI_FLT <- resBHI[which(abs(resBHI[, "log2FoldChange"]) >2& resBHI[,"padj"] < 0.05),]


resMg <- results(dds, contrast=c("group", "coculturewith_Mg", "monoculturewith_Mg"))
resMg_FLT <- resMg[which(abs(resMg[, "log2FoldChange"]) >2& resMg[,"padj"] < 0.05),]

#res_coculture_Mgdependent_FLT <- resBHI[which(abs(resBHI[, "log2FoldChange"]) >1& resBHI[,"padj"] < 0.05& abs(resMg[, "log2FoldChange"]) <2) , ]

res_coculture_Mgdependent_FLT <-  read.csv("DEGs_coculture_Mgdependent.csv")

#candidates_Tnseq <- read.csv("Candidates_Tnseq_BHIY1BHI.csv", header =FALSE, row.names = 1)

candidates_Tnseq <- read.csv("Candidates_Tnseq_BHIY1BHI.csv", header =TRUE)


### make heatmap for interesting genes

heatmap.2(assay(rld)[candidates_Tnseq[[1]],],
          scale="row", trace="none",
          dendrogram="row",
          Colv = NA,
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

#specify clustering and distance function(this setting generates a clear pattern of clustering)
# Reorder the columns of samples based on manuel specification
new_order <- c("BHI.1", "BHI.2", "BHI.3", "BHIY1.1", "BHIY1.2", "BHIY1.3", "BHI.Mg.1", "BHI.Mg.2", "BHI.Mg.3","BHIY1.Mg.1", "BHIY1.Mg.2", "BHIY1.Mg.3")

rld_subset <- assay(rld)[candidates_Tnseq[[1]],]
rld_reordered <- rld_subset[, new_order]


heatmap.2(rld_reordered,
          scale="row", trace="none",
          dendrogram="row",
          Colv = NA,
          #Rowv= NA,
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"),
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))


library("pheatmap")
library("RColorBrewer")

#add group annotation
anno <- data.frame(row.names= candidates_Tnseq$gene, group =candidates_Tnseq$fitness_effect)

#creat color for each group
newCols <- colorRampPalette(grDevices::rainbow(length(unique(anno$group))))
annoCol <- newCols(length(unique(anno$group)))
names(annoCol) <- unique(anno$group)
annoCol <- list(category = annoCol)

pheatmap(rld_reordered, cluster_rows=TRUE, show_rownames=TRUE,scale="row",
         cluster_cols=FALSE, 
         border_color = "black",
         legend_breaks = c(-2, 0, 2),
         clustering_method = "average",
         #cutree_cols = 2,
         color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
         annotation_row = anno, 
         annotation_colors = annoCol)

pheatmap(rld_reordered, cluster_rows=TRUE, show_rownames=TRUE,scale="row",
         cluster_cols=FALSE, 
         border_color = "black",
         legend_breaks = c(-2, 0, 2),
         clustering_method = "average",
         #cutree_cols = 2,
         color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

heatmap.2(assay(rld)[res_coculture_Mgdependent_FLT[[1]],],
          scale="row", trace="none",
          dendrogram="row",
          Colv = NA,
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

heatmap.2(assay(rld)[rownames(resMg_FLT),],
          scale="row", trace="none",
          dendrogram="both" , 
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

## 20240514 wrtie csv file for the raw data counts of Fif S8B.
write.csv (rld_reordered, file ="raw_data_heatmap_FigS8.csv")

###Make a heatpmap and log2FC panel for genes of interest (DEGs between two conditions)
resBHI <- results(dds, contrast=c("group", "cocultureno_Mg", "monocultureno_Mg"))
resBHI_df <- as.data.frame(resBHI)

#Get logsFC of genes of interest
candidates_Tnseq <- read.csv("Candidates_Tnseq_BHIY1BHI.csv", header =TRUE)
resBHI_subset <- resBHI_df[candidates_Tnseq[[1]], ]

Mg_transporter <- read.csv("Mgtransporter.csv", header =TRUE)

resBHI_subset <- resBHI_df[Mg_transporter[[1]], ]

foldchange <- data.frame(row.names = rownames(resBHI_subset), log2FoldChange = resBHI_subset$log2FoldChange)

mean <- data.frame(row.names = rownames(resBHI_subset), mean = resBHI_subset$baseMean)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap", force =TRUE)



library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(foldchange),0, max(foldchange)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean$mean)[1], quantile(mean$mean)[4]), c("white", "red"))

 


rld <- rlog(dds, blind=FALSE)

new_order <- c("BHI.1", "BHI.2", "BHI.3", "BHIY1.1", "BHIY1.2", "BHIY1.3", "BHI.Mg.1", "BHI.Mg.2", "BHI.Mg.3","BHIY1.Mg.1", "BHIY1.Mg.2", "BHIY1.Mg.3")

rld_subset <- assay(rld)[Mg_transporter[[1]],]
rld_reordered <- rld_subset[, new_order]

##create a variable for scaling by row. But re-scaling does not work well. So I give up using this type of heatmap
#scaled_rld = t(scale(t(rld_reordered)))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(rld_reordered,cluster_rows = T, 
              column_labels = colnames(rld_reordered), 
              name="Z-score",
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "average",
              cluster_columns = F)

h2 <- Heatmap(as.matrix(foldchange), row_labels = rownames(rld_reordered), 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(foldchange[i, j],2), x, y)
              })
h3 <- Heatmap(as.matrix(mean), row_labels = rownames(rld_reordered), 
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h

png("./heatmap_v1.png", res = 300, width = 3000, height = 5500)
print(h)
dev.off()
print(h)


## Make a combined heatmap with Tn-seqFC and RNAseq-FC (the below lines don't work well. need to figure them out)
#RNAseq-FC of candidate genes
candidates_Tnseq <- read.csv("Candidates_Tnseq_BHIY1BHI.csv", header =FALSE)
resBHI_subset <- resBHI_df[candidates_Tnseq[[1]], ]

foldchange <- data.frame(row.names = rownames(resBHI_subset), log2FoldChange = resBHI_subset$log2FoldChange)

col_logFC <- colorRamp2(c(min(foldchange),0, max(foldchange)), c("blue", "white", "red"))

#Tnseq FC of candidate genes
Tnseq_foldchange <- read.csv("Tnseq_candidates_FC.csv", header =TRUE)

col_TnseqlogFC <- colorRamp2(c(min(Tnseq_foldchange$log2FoldChange),0, max(Tnseq_foldchange$log2FoldChange)), c("blue", "white", "red")) 

#make a combined heatmap
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h4 <- Heatmap(as.matrix(foldchange), 
              row_labels = rownames(Tnseq_foldchange), 
              cluster_rows = F, name="RNAseq_logFC", 
              top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(foldchange[i, j],2), x, y)
              })

h5 <- Heatmap(as.matrix(Tnseq_foldchange), 
              row_labels = rownames(Tnseq_foldchange), 
              cluster_rows = F, name="Tnseq_logFC", 
              top_annotation = ha, col = col_TnseqlogFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(foldchange[i, j],2), x, y)
              })
h<- h4+h5
h
# heatmap.2(assay(rld_twoThingsWithInteraction)[rownames(res_twoThingsWithInteraction_cult_interestingGenes),],
#           scale="row", trace="none",
#           dendrogram="both" , 
#           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
# 
# heatmap.2(assay(rld_twoThingsWithInteraction)[rownames(res_twoThingsWithInteraction_Mg_interestingGenes),],
#           scale="row", trace="none",
#           dendrogram="both" , 
#           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
# 
# 
# heatmap.2(assay(rld_twoThingsWithInteraction)[rownames(res_twoThingsWithInteraction_interaction_interestingGenes),],
#           scale="row", trace="none",
#           dendrogram="both" , 
#           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
# 
# 
# heatmap.2(assay(rld_twoThingsWithInteraction)[res_twoThingsWithInteraction_interestingGenes_allComparison,],
#           scale="row", trace="none",
#           dendrogram="row" , 
#           col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))







## plot counts of a specific gene across conditions
plotCounts(dds, gene=which.min(resBHI$padj), intgroup="group")

plotCounts(dds, gene="PA4823", intgroup="group")

## Make a volcano plot between BHIY1 and BHI. The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6
library(EnhancedVolcano)

EnhancedVolcano(resMg,
                lab = rownames(resMg),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab =c('PA4824','PA4825','PA4826', 'PA4822', 'PA4823'),
                #ylim =c(0,20), xlim = c(-10,25),
                #pCutoff = 0.05, 
                pCutoffCol = 'padj',
                pointSize= 2,
                subtitle = paste0('p-value cutoff (black line) drawn is 10e-6'
                ),
                title ='MgY1 versus Mg')


###The following script is to make heatmap using InteractiveComplexHeatmap package
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~culture_condition + Mg_condition) 
keep = rowSums(counts(dds)) >= 120
dds = dds[keep, ]

dds$group <- factor(paste0(dds$culture_condition, dds$Mg_condition))

design(dds) <- ~ group
dds$group <- relevel(dds$group, ref ="monocultureno_Mg")
dds <- DESeq(dds)
resultsNames(dds)
interactivate(dds)
res <- results(dds)
res <- as.data.frame (res)
# 
# dds$treatment <- relevel(dds$treatment, ref ="BHI")
# dds <- DESeq(dds)
# res <- results(dds)
# res <- as.data.frame(res)
# 
# res1 <- results(dds)

##make an complex interactive heatmap (from https://github.com/jokergoo/InteractiveComplexHeatmap) I can replace the following lines with interactivate(dds)


#20230510, Janet suggests to put res results with different designs as inputs for plotting heatmap
#res <- as.data.frame(res_twoThingsWithInteraction_cult)
#res <- as.data.frame(res_twoThingsWithInteraction_Mg)
#res <- as.data.frame(res_twoThingsWithInteraction_interaction)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
BiocManager::install("InteractiveComplexHeatmap")

library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)

env = new.env()

make_heatmap = function(fdr = 0.01, base_mean = 0, log2fc = 0) {
  l = res$padj <= fdr & res$baseMean >= base_mean & 
    abs(res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE
    
    if(sum(l) == 0) return(NULL)
    
    m = counts(dds, normalized = TRUE)
    m = m[l, ]
    
    env$row_index = which(l)
    
    ht = Heatmap(t(scale(t(m))), name = "z-score",
                 top_annotation = HeatmapAnnotation(
                   dex = colData(dds)$dex,
                   sizeFactor = anno_points(colData(dds)$sizeFactor)
                 ),
                 show_row_names = FALSE, show_column_names = FALSE, row_km = 2,
                 column_title = paste0(sum(l), " significant genes with FDR < ", fdr),
                 show_row_dend = FALSE) + 
      Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
              name = "log10(baseMean+1)", show_column_names = FALSE) +
      Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
              name = "log2FoldChange", show_column_names = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")))
    ht = draw(ht, merge_legend = TRUE)
    ht
}

# make the MA-plot with some genes highlighted
make_maplot = function(res, highlight = NULL) {
  col = rep("#00000020", nrow(res))
  cex = rep(0.5, nrow(res))
  names(col) = rownames(res)
  names(cex) = rownames(res)
  if(!is.null(highlight)) {
    col[highlight] = "red"
      cex[highlight] = 1
  }
  x = res$baseMean
  y = res$log2FoldChange
  y[y > 2] = 2
  y[y < -2] = -2
  col[col == "red" & y < 0] = "darkgreen"
    par(mar = c(4, 4, 1, 1))
    
    suppressWarnings(
      plot(x, y, col = col, 
           pch = ifelse(res$log2FoldChange > 2 | res$log2FoldChange < -2, 1, 16), 
           cex = cex, log = "x",
           xlab = "baseMean", ylab = "log2 fold change")
    )
}

# make the volcano plot with some genes highlited
make_volcano = function(res, highlight = NULL) {
  col = rep("#00000020", nrow(res))
  cex = rep(0.5, nrow(res))
  names(col) = rownames(res)
  names(cex) = rownames(res)
  if(!is.null(highlight)) {
    col[highlight] = "red"
      cex[highlight] = 1
  }
  x = res$log2FoldChange
  y = -log10(res$padj)
  col[col == "red" & x < 0] = "darkgreen"
    par(mar = c(4, 4, 1, 1))
    
    suppressWarnings(
      plot(x, y, col = col, 
           pch = 16, 
           cex = cex,
           xlab = "log2 fold change", ylab = "-log10(FDR)")
    )
}

## Use shinydashboard to organize the interactive heatmap components
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("shinydashboard")
BiocManager::install("DT")

library(shiny)
library(shinydashboard)
library(DT)

body = dashboardBody(
  fluidRow(
    column(width = 4,
           box(title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               originalHeatmapOutput("ht", height = 800, containment = TRUE)
           )
    ),
    column(width = 4,
           id = "column2",
           box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               subHeatmapOutput("ht", title = NULL, containment = TRUE)
           ),
           box(title = "Output", width = NULL, solidHeader = TRUE, status = "primary",
               HeatmapInfoOutput("ht", title = NULL)
           ),
           box(title = "Note", width = NULL, solidHeader = TRUE, status = "primary",
               htmlOutput("note")
           ),
    ),
    column(width = 4,
           box(title = "MA-plot", width = NULL, solidHeader = TRUE, status = "primary",
               plotOutput("ma_plot")
           ),
           box(title = "Volcanno plot", width = NULL, solidHeader = TRUE, status = "primary",
               plotOutput("volcanno_plot")
           ),
           box(title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
               DTOutput("res_table")
           )
    ),
    tags$style("
            .content-wrapper, .right-side {
                overflow-x: auto;
            }
            .content {
                min-width:1500px;
            }
        ")
  )
)

## brush action on a heatmap
library(DT)
library(GetoptLong) # for qq() function
brush_action = function(df, input, output, session) {
  
  row_index = unique(unlist(df$row_index))
  selected = env$row_index[row_index]
  
  output[["ma_plot"]] = renderPlot({
    make_maplot(res, selected)
  })
  
  output[["volcanno_plot"]] = renderPlot({
    make_volcano(res, selected)
  })
  
  output[["res_table"]] = renderDT(
    formatRound(datatable(res[selected, c("baseMean", "log2FoldChange", "padj")], rownames = TRUE), columns = 1:3, digits = 3)
  )
  
  output[["note"]] = renderUI({
    if(!is.null(df)) {
      HTML(qq("<p>Row indices captured in <b>Output</b> only correspond to the matrix of the differential genes. To get the row indices in the original matrix, you need to perform:</p>
<pre>
l = res$padj <= @{input$fdr} & 
    res$baseMean >= @{input$base_mean} & 
    abs(res$log2FoldChange) >= @{input$log2fc}
l[is.na(l)] = FALSE
which(l)[row_index]
</pre>
<p>where <code>res</code> is the complete data frame from DESeq2 analysis and <code>row_index</code> is the <code>row_index</code> column captured from the code in <b>Output</b>.</p>"))
    }
  })
}

## Side bar information
ui = dashboardPage(
  dashboardHeader(title = "DESeq2 results"),
  dashboardSidebar(
    selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05)),
    numericInput("base_mean", label = "Minimal base mean:", value = 0),
    numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0),
    actionButton("filter", label = "Generate heatmap")
  ),
  body
)

## On the server side
server = function(input, output, session) {
  observeEvent(input$filter, {
    ht = make_heatmap(fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc)
    if(!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
                                    brush_action = brush_action)
    } else {
      # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap, thus, it is ht_heatmap here.
      output$ht_heatmap = renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui, server)


