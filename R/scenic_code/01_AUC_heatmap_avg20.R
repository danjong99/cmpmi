library(corrplot)
library(gplots)
library(tidyr)
library(dplyr)
library(NMF)
library(RColorBrewer)
library(reshape2)

### Extract AUC matrix
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
auc_mtx <- getAUC(aucell_regulonAUC)
saveRDS(auc_mtx, file = "auc_mtx.Rds")

### Load from the saved RDS
auc_mtx <- readRDS("./R/scenic_code/rdata/auc_mtx.Rds")
### Load Seurat RDS
cmp <- readRDS("./RDS/seurat_combined.rds")
cmp <- UpdateSeuratObject(cmp)
cell_index <- cmp@meta.data
cell_index$res.1 <- NULL
cell_index$cellTrack <- rownames(cell_index)
cell_index$cellTrack <- gsub("WT_","",cell_index$cellTrack)
cell_index$cellTrack <- gsub("GF_","",cell_index$cellTrack)
sum(duplicated(cell_index$cellTrack))
sum(duplicated(intersect(colnames(auc_mtx), cell_index$cellTrack)))
cell_index[cell_index$cellTrack == colnames(auc_mtx)[1],]

a <- list()
for (i in 1:length(colnames(auc_mtx))) {
  a[i] <- grep(colnames(auc_mtx)[i], cell_index$cellTrack)
}

avg20_cmp <- subset(cmp, cells = rownames(cell_index[unlist(a),]))
colnames(auc_mtx) <- rownames(cell_index[unlist(a),])
avg20_cmp <- avg20_cmp[,colnames(auc_mtx)]
cellInfo <- avg20_cmp@meta.data
cellInfo$nGene <- NULL; cellInfo$nUMI <- NULL; cellInfo$res.1 <- NULL
cellInfo$orig.celltype <- cellInfo$res.1.2
colnames(cellInfo) <- c("orig.ident","CellType","orig.celltype")
cellInfo$CellType <- gsub("0","00",cellInfo$CellType)
cellInfo$CellType <- gsub("1","01",cellInfo$CellType)
cellInfo$CellType <- gsub("2","02",cellInfo$CellType)
cellInfo$CellType <- gsub("3","03",cellInfo$CellType)
cellInfo$CellType <- gsub("4","04",cellInfo$CellType)
cellInfo$CellType <- gsub("5","05",cellInfo$CellType)
cellInfo$CellType <- gsub("6","06",cellInfo$CellType)
cellInfo$CellType <- gsub("7","07",cellInfo$CellType)
cellInfo$CellType <- gsub("8","08",cellInfo$CellType)
cellInfo$CellType <- gsub("9","09",cellInfo$CellType)
cellInfo$CellType <- gsub("0100","10",cellInfo$CellType)
cellInfo$CellType <- gsub("0101","11",cellInfo$CellType)
cellInfo$CellType <- gsub("0102","12",cellInfo$CellType)
cellInfo <- cellInfo[order(cellInfo$CellType, decreasing = FALSE),]
cellInfomac <- subset(cellInfo, cellInfo$CellType == "01"|cellInfo$CellType == "02"|cellInfo$CellType == "03"|
                        cellInfo$CellType == "04"|cellInfo$CellType == "06"|cellInfo$CellType == "07"|cellInfo$CellType == "11")
auc_mtx <- as.data.frame(auc_mtx)[rownames(cellInfo)]

#Heatmap
a <- auc_mtx
cData <- t(apply(a,1,function(x)x-median(x)))
scale.range <- c(-2,2)
scale.brakes <- seq(scale.range[1],scale.range[2], by = 0.05)
n.brakes <- length(scale.brakes) - 1
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(n.brakes)
colVars <- readRDS("./rdata/colVars.Rds")

rowDendro <- as.dendrogram(hclust(as.dist(1-cor(t(cData)))))
colDendro <- as.dendrogram(hclust(as.dist(1-cor(cData))))
heatmap.2(cData, trace = 'none', scale = 'row', labRow = rownames(cData), labCol = F,
          col = Colors, breaks = scale.brakes, symkey = T, keysize = 1,
          cexRow = 0.3, Rowv = rowDendro, Colv = colDendro, dendrogram = c("both"))

png(filename = "GENIE3_TFs_Avg20_4.tiff",
    width = 8000, height = 6000, units = "px", pointsize = 9, res = 600)
heatmap.2(cData, trace = 'none', scale = 'row', labRow = rownames(cData), labCol = F,
          col = Colors, breaks = scale.brakes, symkey = T, keysize = 1,
          cexRow = 0.3,
          Rowv = rowDendro, Colv = "NA", dendrogram = c("row"),
          ColSideColors = c(
            rep(colVars[["CellType"]][["0"]],sum(cellInfo$CellType == "00")),
            rep(colVars[["CellType"]][["1"]],sum(cellInfo$CellType == "01")),
            rep(colVars[["CellType"]][["2"]],sum(cellInfo$CellType == "02")),
            rep(colVars[["CellType"]][["3"]],sum(cellInfo$CellType == "03")),
            rep(colVars[["CellType"]][["4"]],sum(cellInfo$CellType == "04")),
            rep(colVars[["CellType"]][["5"]],sum(cellInfo$CellType == "05")),
            rep(colVars[["CellType"]][["6"]],sum(cellInfo$CellType == "06")),
            rep(colVars[["CellType"]][["7"]],sum(cellInfo$CellType == "07")),
            rep(colVars[["CellType"]][["8"]],sum(cellInfo$CellType == "08")),
            rep(colVars[["CellType"]][["9"]],sum(cellInfo$CellType == "09")),
            rep(colVars[["CellType"]][["10"]],sum(cellInfo$CellType == "10")),
            rep(colVars[["CellType"]][["11"]],sum(cellInfo$CellType == "11")),
            rep(colVars[["CellType"]][["12"]],sum(cellInfo$CellType == "12")))
)
dev.off()

#Add CellType to Melt Matrix
auc_mtx$Regulons <- rownames(auc_mtx)
auc_mtx_melt <- melt(auc_mtx, value.name = "Regulons")
auc_mtx_melt$CellType <- "00"
colnames(auc_mtx_melt) <- c("Regulons","cells","AUC","CellType")
for (i in 1:length(auc_mtx_melt$cells)) {
  auc_mtx_melt$CellType[i] <- cellInfo[rownames(cellInfo) == auc_mtx_melt$cells[i],]$CellType
}
head(auc_mtx_melt)
tail(auc_mtx_melt)

#Average AUC Valuse of each regulon from every celltype
auc_mtx_melt_summarize <- auc_mtx_melt %>% group_by(Regulons, CellType) %>% summarize(avg_AUC = mean(AUC))
head(auc_mtx_melt_summarize)
auc_mtx_average <- as.matrix(acast(auc_mtx_melt_summarize, Regulons ~ CellType, value.var="avg_AUC", fill="none"))
head(auc_mtx_average)

# TF correlation
class(auc_mtx_average) <- "numeric"
auc_mtx_average <- round(auc_mtx_average, 3)
M <- cor((auc_mtx_average))
res1 <- cor.mtest((auc_mtx_average), conf.level = 0.95)
res2 <- cor.mtest((auc_mtx_average), conf.level = 0.99)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
wb <- c("white", "black")

#### Closeness of clusters based on Averaged AUC
corrplot(M, method="color", order = "hclust", hclust.method = "complete",
         tl.col = "black", col = col2(100), tl.cex = 0.5, 
         addrect = 3, rect.col = "red", rect.lwd = 3, 
         addshade = c("negative"),shade.col = "#053061", cl.pos = "b")