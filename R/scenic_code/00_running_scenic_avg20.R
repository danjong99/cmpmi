######## SCENIC should be run on Server Computer

library(Seurat)
library(SCENIC)
library(SingleCellExperiment)
library(dplyr)

set.seed(124)
avg20 <- read.csv("./source/avg20.csv")

rownames(avg20) <- avg20$GENE
avg20$GENE <- NULL
avg20 <- as.matrix(avg20)
avg20 <- avg20*20
head(avg20)

cmp <- readRDS("seurat_combined.rds")
cellLabels <- cmp@meta.data
rm(cmp)
cellLabels$cellID <- rownames(cellLabels)
cellLabels <- cellLabels[,c("res.1.2","cellID")]
cellLabels$cellID <- gsub("WT_","",cellLabels$cellID)
cellLabels$cellID <- gsub("GF_","",cellLabels$cellID)
cellLabels <- cellLabels[!duplicated(cellLabels$cellID),]
rownames(cellLabels) <- cellLabels$cellID
colnames(cellLabels) <- c("CellType","cellID")

index <- intersect(rownames(cellLabels),colnames(avg20))
cellLabels <- cellLabels[index,]
cellLabels$cellID <- NULL
avg20 <- avg20[,index]
length(colnames(avg20))

avg20_sceMP <- SingleCellExperiment(assays = list(counts = avg20),
                                     colData=data.frame(cellLabels[colnames(avg20),, drop=FALSE]))
dir.create("rdata")
save(avg20_sceMP, file="./rdata/avg20_sceMP.RData")

##### you can start from here
load("./rdata/avg20_sceMP.RData")
exprMat <- counts(avg20_sceMP)
dim(exprMat)

cellInfo <- colData(avg20_sceMP)
cellInfo$nGene <- colSums(exprMat>0)
cellInfo <- data.frame(cellInfo)
head(cellInfo)
setwd("./")
dir.create("int")
saveRDS(cellInfo, file="./int/cellInfo.Rds")
org="mgi" # or hgnc, or dmel
dbDir="path_to_databases" # RcisTarget databases location
myDatasetTitle="scenic_colon_mac" # choose a name for your analysis
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=4) 
scenicOptions@inputDatasetInfo$cellInfo <- "./int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "./int/colVars.Rds"
saveRDS(scenicOptions, file="./int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))
interestingGenes <- c("Runx3", "Atf3", "Hes1")
interestingGenes[which(!interestingGenes %in% genesKept)]
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

corrMat <- cor(t(exprMat_filtered), method="spearman")
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
set.seed(124)
exprMat_filtered <- log2(exprMat_filtered+1)

######## Wrapper Functions of Scenic ###################
runGenie3(exprMat_filtered, scenicOptions)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 124

runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions@settings[["dbDir"]] <- dbDir
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,15,50), perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,15,50), perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/): 
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE)

par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
plotTsne_compareSettings("int/tSNE_AUC_50pcs_50perpl.Rds", scenicOptions, showLegend=FALSE, varName="CellType", cex=2)

par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 50
scenicOptions@settings$defaultTsne$perpl <- 50
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

runSCENIC_4_aucell_binarize(scenicOptions)
########################################################################

logMat <- exprMat_filtered
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) #default t-SNE
savedSelections <- shiny::runApp(aucellApp)

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
library(AUCell)
auc_mtx <- getAUC(aucell_regulonAUC)

# Show TF expression:
dev.off()
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Elf2")],], plots="Expression")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Elf2")],], plots="AUC")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Crem")],], plots="Expression")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Crem")],], plots="AUC")

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, exprMat = exprMat_filtered, plots="Expression")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots = "AUC")

library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

par(mfrow=c(1,2))
regulonNames <- c("Atf3")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
regulonNames <- list(red=c("Atf3"),
                     green=c("Klf5"),
                     blue=c( "Irf7"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(-30,-25-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
regulonTargetsInfo <- RcisTarget::addLogo(regulonTargetsInfo, motifCol="bestMotif")
regulonTargetsInfo$Genie3Weight <- signif(regulonTargetsInfo$Genie3Weight, 2)
View(regulonTargetsInfo[regulonTargetsInfo$gene == "Tubb5",])

library(DT)
colsToShow <- c("TF", "gene", "nMotifs", "bestMotif", "logo", "NES", "highConfAnnot", "Genie3Weight")
DT::datatable(regulonTargetsInfo[TF=="" & highConfAnnot==TRUE, colsToShow, with=F], escape=FALSE, filter="top")
DT::datatable(regulonTargetsInfo[highConfAnnot==TRUE, colsToShow, with=F], escape=FALSE, filter="top")
DT::datatable(regulonTargetsInfo[colsToShow], escape=FALSE, filter="top")

TFs_lookat <- c("Crem","Prdm1","Junb","Cebpb","Atf4","Jund","Klf9","Egr1","Jun","Egr2","Irf1")
tfs_selected <- data.frame()
b <- data.frame()
for (i in 1:length(TFs_lookat)) {
  b <- tf_target_table %>% filter(TF == TFs_lookat[i]) %>% select("TF","gene","bestMotif","Genie3Weight")
  tfs_selected <- rbind(tfs_selected, b)
}
tfs_selected <- tfs_selected %>% filter(TF != gene)
write.csv(tfs_selected, "tf_targets_selected.csv", row.names = FALSE)

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
motifEnrichment_selfMotifs_wGenes <- RcisTarget::addLogo(motifEnrichment_selfMotifs_wGenes)
colsToShow <- c("motifDb", "logo", "NES", "geneSet", "TF_highConf") # "TF_lowConf", "enrichedGenes"
DT::datatable(motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Pbx1", colsToShow, with=F], escape=FALSE, filter="top")

getOutName(scenicOptions, "s2_regulonTargetsInfo")

avg20.4.index <- read.delim("~/Desktop/10x_scRNAseq/GF1_WT2_MergedRCC/SCENIC_R/Avg20_All/No4/source/avg20.4.index")
avg20.4.index <- t(avg20.4.index)
avg20.4.index <- cbind(rownames(avg20.4.index), avg20.4.index)
rownames(avg20.4.index) <- NULL

