## This is the script I use for analyzing data of BHI_Tnseq_20221020. This script is modified from Tnseq_code_for_Phoebe_2021_Oct
library(dplyr) ## to manipulate tables 
library(GenomicRanges)  ## to look at gene overlaps

#Clear environment
rm(list=ls())

## define some files and locations
# my mac does not mount the shared drive in a consistent location - annoying! 
# so I have first figure out where it is:
sharedDriveLocation <- "/Volumes/fh/fast/malik_h/"
if (!file.exists(sharedDriveLocation)) { sharedDriveLocation <- "/Volumes/malik_h/" }

## change working dir
workingDir <- paste(sharedDriveLocation, "/user/yhsieh/Raw_Sequence_Reads/Tnseq_20221020/usftp21.novogene.com/R_analysis", sep="")
setwd(workingDir)

## specify where Phoebe's analysis files are 
inputFileDir <- paste(sharedDriveLocation, "", sep="/user/yhsieh/Raw_Sequence_Reads/Tnseq_20221020/usftp21.novogene.com/R_analysis")

## get filenames of the sites files for each sample (*-sites.txt). These sites files are the output of one of the Whiteley lab's scripts. 
# two columns - number of reads at that site, position in genome. Sorted with most common sites at the top
# They are space delimited, also have spaces at the start of each line
# I edit the following line to read sites.txt files in the folder of R_analysis in the above directory
sitesFiles <- list.files(inputFileDir, pattern="sites.txt$", full.names = TRUE)

## make a custom R function to read in a sites file
# we'll add column names and provide an option to remove the most common sites in each sample
myReadSitesFunction <- function(oneSitesFile, removeMostCommonSites=NULL) {
  ## read in the file. The read.table function can deal with the space-delimited file
  sites <- read.table(oneSitesFile)
  
  ## add column names.  For the column showing number of reads, we will name the column using the sample name, which we get from the file name
  sampleName <- strsplit(oneSitesFile, "/")[[1]] # split by "/"
  sampleName <- sampleName[length(sampleName)]   # take the last item
  sampleName <- gsub("-sites.txt","",sampleName)
  colnames(sites) <- c(sampleName, "position")
  
  ## if removeMostCommonSites is not null, then we remove the most common sites
  # the sites table is ALREADY sorted with the most common sites at the top, so we will remove the first however many rows
  if (!is.null(removeMostCommonSites)) {
    startingRow <- removeMostCommonSites+1
    endingRow <- dim(sites)[1]
    sites <- sites[startingRow:endingRow,]
  }
  ## we're done - return the output
  return(sites)
}

## if we wanted to run that function on one file, we would do this:
# temp <- myReadSitesFunction(sitesFiles[1], removeMostCommonSites=50)

## run that function on all 7 sites files, removing the 50 most common sites in each sample
# lapply is a VERY useful function!  it is worth learning about it. lapply is a way to apply a customized function to a list
# this step is to remove the top 50 sites with the most reads for normalizing amplification bias
# the result is a list object, each item is a data.frame
sitesTables <- lapply(sitesFiles, myReadSitesFunction, removeMostCommonSites=50)

## get an object that's just the sample names - it will be useful later.  We get it from the first column name of each table
# the sapply function is very similar to lapply, but instead of returning a list object as a result, it tries to simplify the output(vector or matrix) - e.g.
sampleNames <- sapply(sitesTables, function(x) { colnames(x)[1] }) 

# also add names to the sitesTables list
names(sitesTables) <- sampleNames

## show how many rows and columns in each table:
# t function transposes a table (swaps rows and columns)
t(sapply(sitesTables, dim))



###20211117-to make sites.txt(counts, position) into an input file(only consider the start position of reads right next to Tn) for IGV. wrtieBEDGRAPHfile converts each start site into Bedgraph format with 1bp from the start(so called "zero-based half-open"). I can load Bedgraph file of each sample to IGV to check Tn insertion in genes of interest. Worth thinking about if I need to have a bedgraph file representing ratio between 2 samples

# filename: output file of BEDGRAPH function
writeBEDGRAPHfile <- function(oneSiteTable, filename, normalize=FALSE) {
  newTable <- data.frame(chr="refseq|NC_002516.2|chromosome",
                         start=(oneSiteTable[,"position"]-1),
                         end=oneSiteTable[,"position"],
                         value=oneSiteTable[,1])
  # maybe normalize
  if(normalize) {
    newTable[,"value"] <- newTable[,"value"] / (sum(newTable[,"value"])/1000000)
  }
  # sort by position
  newTable <- newTable[ order(newTable[,"start"]) , ]
  # write the file
  write.table(newTable, file=filename, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}

for (thisSample in names(sitesTables)) {
  thisFileName <- paste("Bedgraph/",thisSample,"_sites.bedGraph",sep="")
  cat("writing file",thisFileName,"\n")
  writeBEDGRAPHfile(sitesTables[[thisSample]], thisFileName)
}


# This section outputs normalized bedgraph files
for (thisSample in names(sitesTables)) {
  thisFileName <- paste("Bedgraph/",thisSample,"_normsites.bedGraph",sep="")
  cat("writing normalized file",thisFileName,"\n")
  writeBEDGRAPHfile(sitesTables[[thisSample]], thisFileName, normalize = TRUE)
}




######### now we want to combine all those sites tables into one giant table, adding 0s in samples where a site was not seen. We will also sort by genomic position. There is probably a more elegant way to do this but this works just fine

## first make a function to change any NA value to 0
replaceNAwith0 <- function(x) {
  x [ which(is.na(x)) ] <- 0
  return(x)
}

### then merge all the tables (all have a column named "position")
# the full_join function figures out columns that are shared across tables and uses those to merge
sitesTables_merged <- sitesTables[[1]][,c(2,1)] %>% 
  full_join(sitesTables[[2]]) %>% 
  full_join(sitesTables[[3]]) %>% 
  full_join(sitesTables[[4]]) %>% 
  full_join(sitesTables[[5]]) %>% 
  full_join(sitesTables[[6]]) %>% 
  full_join(sitesTables[[7]]) %>% 
  full_join(sitesTables[[8]]) %>%
  full_join(sitesTables[[9]]) %>%
  full_join(sitesTables[[10]]) %>% 
  full_join(sitesTables[[11]]) %>% 
  full_join(sitesTables[[12]]) %>% 
  full_join(sitesTables[[13]]) %>% 
  full_join(sitesTables[[14]]) %>%
  full_join(sitesTables[[15]]) %>% 
  full_join(sitesTables[[16]]) %>% 
  full_join(sitesTables[[17]]) %>% 
  arrange(position) %>%  # sort the resulting merged table by position
  mutate(across(all_of(sampleNames), replaceNAwith0) ) # for each sample column, replace NAs with 0

## add a column that says how many samples each site was seen in:
# the apply function is also useful. Works on a data.frame to apply a function to every row or column
sitesTables_merged <- sitesTables_merged %>% 
  mutate( numSamplesWithSite = 
            apply(sitesTables_merged[,sampleNames], 1, function(x) {
              sum(x>0)
            }))

### 20211117-the following is to make an input file for DEseq2, using sites table. I rename the rows to be the start position of each read. This will allow me to do comparison and statistical analysis per sites. After getting the list of candidate sites, I can match information of gene vs CDS with the sites.
rownames(sitesTables_merged) <-paste (sitesTables_merged[,"position"], "bp", sep="")
save(sitesTables_merged, file ="sitesTables_merged.rda")

###### because we will use some functions from the GenomicRanges Bioconductor package to understand how genes and sites overlap, we will make an 'IRanges' style-object that corresponds to sitesTables_merged:
sitesTables_merged_ir <- IRanges(start=sitesTables_merged[,"position"], width=1)
values(sitesTables_merged_ir) <- sitesTables_merged[,sampleNames]

head(sitesTables_merged_ir)

##### read in gene coordinates from the gff file
geneCoordsFile <- paste(sharedDriveLocation, "user/yhsieh/PAO1_Reference_Genome/Pseudomonas_aeruginosa_PAO1_107_truncated.gff.txt", sep="")

# got this code from the Whiteley script
geneCoords <- read.delim(file=geneCoordsFile, 
                         sep="\t", fill=TRUE, header=FALSE, 
                         col.names = c("seqname", "source", "feature", 
                                       "start", "end", "score", "strand", "frame", "att"))
geneCoords <- geneCoords[(geneCoords$feature=="gene"),]


### I also want to parse the attributes column into separate columns so they're easier to look at
temp <- strsplit(geneCoords[,"att"], ";")
# each gene has either 3 or 4 attributes
table(sapply(temp, length))
#    3    4 
# 3430 2278 

# get the first, second and third fields and make new columns
geneCoords[,"ID"] <- gsub("ID=","",sapply(temp, "[[", 1))
geneCoords[,"Alias"] <- gsub("Alias=","",sapply(temp, "[[", 2))
geneCoords[,"Dbxref"] <- gsub("Dbxref=PGD_Gene_ID:","",sapply(temp, "[[", 3))

# for the fourth field (name) we do something more complicated because not all genes have it
geneCoords[,"name"] <- NA
tempIndicesOfGenesWithNameField <- which(sapply(temp, length)==4)
geneCoords[tempIndicesOfGenesWithNameField,"name"] <- 
  gsub("name=","",sapply( temp[tempIndicesOfGenesWithNameField], "[[", 4))

rm(temp, tempIndicesOfGenesWithNameField) # tidy up

### it turns out the Dbxref field is identical to ID so we can get rid of it. We can also get rid of a few other columns that are not very useful, including att, because we parsed all that information and put it in other columns
# check those columns really are identical
table(geneCoords[,"ID"] == geneCoords[,"Dbxref"], useNA="always")
# TRUE <NA> 
# 5708    0 

geneCoords <- geneCoords %>% 
  select(-seqname, -source, -feature, -score, -frame, -att, -Dbxref) # get rid of a few columns 


## like we did for sites_merged, we make a Bioconductor-compatible IRanges ('ir') object from geneCoords
geneCoords_ir <- IRanges(start=geneCoords[,"start"], end=geneCoords[,"end"])
values(geneCoords_ir) <- geneCoords[,c("strand", "ID", "Alias", "name")]

# the following script is added to extract Tn site information of PA0763 
names(geneCoords_ir) <- geneCoords[,"Alias"]
geneCoords_ir['PA0763']
subsetByOverlaps(sitesTables_merged_ir, geneCoords_ir['PA0763'])

#head(geneCoords_ir)
# geneCoords_ir is a very beautiful table


######### get site and read counts per gene for each sample

## mergeByOverlaps is a function from the GenomicRanges packages, that uses the coordinates to figure out which genes and sites overlap each other, and merges our two tables based on the overlaps
overlapTable <- mergeByOverlaps(geneCoords_ir, sitesTables_merged_ir)
head(overlapTable, 15)

## the ID column of overlapTable is a plain character vector, but we'll turn it into a factor. We do this so that R knows which genes are in the genome but not present in the overlaps table, so that we also get counts for genes with 0 insertions observed
overlapTable[,"ID"] <- factor(overlapTable[,"ID"], levels=geneCoords[,"ID"])


## now we're going to summarize those counts for each gene. Start with a table that's just the gene info. 
geneCoordsPlusCounts <- geneCoords 

## add a column to show how many insertion sites were observed in any sample in each gene
# this table has one row per insertion, so we can simply count how many times each gene name appears to count how many insertions it has:
temp <- table(overlapTable[,"ID"])
head(temp)
# check the counts are in the same order as the gene table - they should be:
table( names(temp) == geneCoordsPlusCounts[,"ID"])
geneCoordsPlusCounts[,"numSitesAnySample"] <- temp
rm(temp)

## for each sample, count number of sites observed in each gene
# for loops are also useful!  lapply is similar and a bit more elegant
for (tempSample in sampleNames) {
  tempNewColumnName <- paste(tempSample, "_numSites",sep="")
  # tapply is also very useful - worth learning about this function!  
  # Here we look at the values in the tempSample column of overlapTable, we group them by gene ID, and for each gene ID we count how many values are >0
  # the dplyr equivalent function of tapply is called 'summarise'
  numSitesPerGene <- tapply( overlapTable[,tempSample], overlapTable[,"ID"], function(x) {
    sum(x>0)
  })
  numSitesPerGene <- replaceNAwith0(numSitesPerGene)
  geneCoordsPlusCounts[,tempNewColumnName] <- numSitesPerGene
  rm(numSitesPerGene, tempNewColumnName)
}

## for each sample, get total number of reads observed in each gene
for (tempSample in sampleNames) {
  tempNewColumnName <- paste(tempSample, "_numReads",sep="")
  numReadsPerGene <- tapply( overlapTable[,tempSample], overlapTable[,"ID"], sum)
  numReadsPerGene <- replaceNAwith0(numReadsPerGene)
  geneCoordsPlusCounts[,tempNewColumnName] <- numReadsPerGene
  rm(numReadsPerGene, tempNewColumnName)
}


## get normalized read counts for each gene in sample, by dividing each by total number of mapped reads.  

# I'll get total number of reads mapped anywhere in the genome for each sample, in millions (actually, this does ignore the 50 most common sites, which we removed before making sitesTables_merged). The Whiteley script does this slightly differently - they use DESeq's EstimateSizeFactors on the per-gene counts, so they are ignoring any intergenic sites in their totals for each sample. Not sure which is better.
normalizationFactors <- apply(sitesTables_merged[,sampleNames], 2, sum) / 1000000

# now divide by the relevant normalizationFactors
for (tempSample in sampleNames) {
  tempOldColumnName <- paste(tempSample, "_numReads",sep="")
  tempNewColumnName <- paste(tempSample, "_numReadsNorm",sep="")
  geneCoordsPlusCounts[,tempNewColumnName] <- 
    (geneCoordsPlusCounts[,tempOldColumnName] / normalizationFactors[tempSample])+1
  rm(tempOldColumnName, tempNewColumnName)
}


## 20221108 plot number of normalized Tn insertions per gene between b86_BHI_tet(8) vs b86_BHI_Y1_tet(9),b86_BHI(11) vs B86_BHI_Y1(10)
library(scales)
plot(geneCoordsPlusCounts[,8], geneCoordsPlusCounts[,9], cex=0.25, log = "xy", xlab="b86_BHI_tet", ylab= "b86_BHI_Y1_tet", col =alpha("black", 0.25))

plot(log(geneCoordsPlusCounts[,47],2), log(geneCoordsPlusCounts[,52],2), cex=0.25,  xlab="Tn12_BHIY1", ylab= "b86_BHI_Y1", col =alpha("black", 0.25))
abline(0,1, col="red")
abline(log(2,2),1, col= "orange")
abline(-log(2,2),1, col= "blue")

# compare b86BHIY1 and Tn12_BHIY1 to find interesting genes. 
geneCoordsPlusCounts[,"b86_BHIY1_Tn12_Log2Ratio"]<-
  log2(geneCoordsPlusCounts[, paste(sampleNames[3], "_numReadsNorm",sep="")]/ geneCoordsPlusCounts[, paste(sampleNames[15], "_numReadsNorm",sep="")] )

geneCoordsPlusCounts[,"b86BHIY1_Tn12_interesting"] <- abs(geneCoordsPlusCounts[,"b86_BHIY1_Tn12_Log2Ratio"]) >=1 & 
  (geneCoordsPlusCounts[, paste(sampleNames[3], "_numReadsNorm",sep="")]+ geneCoordsPlusCounts[,paste(sampleNames[15], "_numReadsNorm",sep="")]) >= 10


InterestingGenes_b86BHIY1_Tn12 <- geneCoordsPlusCounts[(geneCoordsPlusCounts[, "b86BHIY1_Tn12_interesting"]),]
write.csv(InterestingGenes_b86BHIY1_Tn12, row.names=FALSE, file= "MappedReads_R/InterestingGenes_b86BHIY1_Tn12.csv")

# compare b86BHIY1 and b86_BHI to find interesting genes. 
geneCoordsPlusCounts[,"b86_BHIY1_b86BHI_Log2Ratio"]<-
  log2(geneCoordsPlusCounts[, paste(sampleNames[3], "_numReadsNorm",sep="")]/ geneCoordsPlusCounts[, paste(sampleNames[4], "_numReadsNorm",sep="")] )

geneCoordsPlusCounts[,"b86BHIY1_b86BHI_interesting"] <- abs(geneCoordsPlusCounts[,"b86_BHIY1_b86BHI_Log2Ratio"]) >=1 & 
  (geneCoordsPlusCounts[, paste(sampleNames[3], "_numReadsNorm",sep="")]+ geneCoordsPlusCounts[,paste(sampleNames[4], "_numReadsNorm",sep="")]) >= 10


InterestingGenes_b86BHIY1_b86BHI <- geneCoordsPlusCounts[(geneCoordsPlusCounts[, "b86BHIY1_b86BHI_interesting"]),]
write.csv(InterestingGenes_b86BHIY1_b86BHI, row.names=FALSE, file= "MappedReads_R/InterestingGenes_b86BHIY1_b86BHI.csv")

# compare b86BHIY1_tet and b86_BHI_tet to find interesting genes. 
geneCoordsPlusCounts[,"b86_tetY1_b86tet_Log2Ratio"]<-
  log2(geneCoordsPlusCounts[, paste(sampleNames[2], "_numReadsNorm",sep="")]/ geneCoordsPlusCounts[, paste(sampleNames[1], "_numReadsNorm",sep="")] )

geneCoordsPlusCounts[,"b86tetY1_b86tet_interesting"] <- abs(geneCoordsPlusCounts[,"b86_tetY1_b86tet_Log2Ratio"]) >=1 & 
  (geneCoordsPlusCounts[, paste(sampleNames[2], "_numReadsNorm",sep="")]+ geneCoordsPlusCounts[,paste(sampleNames[1], "_numReadsNorm",sep="")]) >= 10


InterestingGenes_b86tetY1_b86tet <- geneCoordsPlusCounts[(geneCoordsPlusCounts[, "b86tetY1_b86tet_interesting"]),]
write.csv(InterestingGenes_b86tetY1_b86tet, row.names=FALSE, file= "MappedReads_R/InterestingGenes_b86tetY1_b86tet.csv")

##20210915  save a giant table for input to DEseq2
save (geneCoordsPlusCounts, file ="geneCoordsPlusCounts.rda" )
