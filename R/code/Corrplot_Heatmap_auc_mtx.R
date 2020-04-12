library(corrplot)
library(gplots)
library(tidyr)
library(dplyr)

# Color Choice
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

wb <- c("white", "black")

#Correlation & Confidence Intervals
M<- cor(t(a))
res1 <- cor.mtest(auc_mtx, conf.level = 0.95)
res2 <- cor.mtest(auc_mtx, conf.level = 0.99)
head(round(M,2))
set.seed(0)

corrplot(M, method="color", order = "hclust", 
         tl.col = "white", col = col1(10),
         addrect = 5, rect.col = "white")

corrplot(M, method="color", order = "hclust", hclust.method = c("ward.D2"),
         tl.col = "white", col = col1(100), tl.cex = 0.5)

cellInfo$CellType <- gsub("1","01",cellInfo$CellType)
cellInfo$CellType <- gsub("2","02",cellInfo$CellType)
cellInfo$CellType <- gsub("3","03",cellInfo$CellType)
cellInfo$CellType <- gsub("4","04",cellInfo$CellType)
cellInfo$CellType <- gsub("6","06",cellInfo$CellType)
cellInfo$CellType <- gsub("7","07",cellInfo$CellType)
cellInfo$CellType <- gsub("0101","11",cellInfo$CellType)

index <- rownames(cellInfo[order(cellInfo$CellType, decreasing = FALSE),])
d <- auc_mtx[,index]
identical(colnames(d), colnames(auc_mtx))
length(intersect(colnames(d), colnames(auc_mtx)))

##Heatmap
a <- auc_mtx
cData <- t(apply(a,1,function(x)x-median(x)))
b <- t(cData)[,!(colSums(t(cData)) == 0)]
cData <- t(b)

#cData <- a
scale.range <- c(-2,2)
scale.brakes <- seq(scale.range[1],scale.range[2], by = 0.5)
n.brakes <- length(scale.brakes) - 1
Palette <- colorpanel(n.brakes,'lightcyan1','lightgoldenrodyellow','indianred2')

rowDendro <- as.dendrogram(hclust(as.dist(1-cor(t(cData)))))
colDendro <- as.dendrogram(hclust(as.dist(1-cor(cData))))

heatmap.2(cData, trace = 'none', scale = 'row', labRow = rownames(cData), labCol = F,
          col = Palette, breaks = scale.brakes, symkey = T, keysize = 1, cexRow = 0.5,
          Rowv = rowDendro, Colv = colDendro, dendrogram = c("both"))

rowDendro <- as.dendrogram(hclust(as.dist(1-cor(t(cData)))))
colDendro <- as.dendrogram(hclust(as.dist(1-cor(cData))))
tiff(filename = "Heatmap_TFs_AUC_GRNBOOST.tiff",
     width = 4500, height = 3500, units = "px", pointsize = 12, res = 300)
heatmap.2(cData, trace = 'none', scale = 'row', labRow = rownames(cData), labCol = F,
          col = Palette, breaks = scale.brakes, symkey = T, keysize = 1, cexRow = 0.5,
          Rowv = rowDendro, Colv = "NA", dendrogram = c("row"),
          ColSideColors = c(
            rep(colVars_mac[["CellType"]][["1"]],sum(cellInfo$CellType == "01")),
            rep(colVars_mac[["CellType"]][["2"]],sum(cellInfo$CellType == "02")),
            rep(colVars_mac[["CellType"]][["3"]],sum(cellInfo$CellType == "03")),
            rep(colVars_mac[["CellType"]][["4"]],sum(cellInfo$CellType == "04")),
            rep(colVars_mac[["CellType"]][["6"]],sum(cellInfo$CellType == "06")),
            rep(colVars_mac[["CellType"]][["7"]],sum(cellInfo$CellType == "07")),
            rep(colVars_mac[["CellType"]][["11"]],sum(cellInfo$CellType == "11")))
            )
dev.off()
#heatmap.2(cData, trace = 'none', scale = 'row', labRow = rownames(a), offsetRow = 0.01, 
#          cexRow = 0.5,col = Palette, #breaks = scale.brakes, 
#          symkey = T,
#          Rowv = FALSE, Colv = "NA", dendrogram = c("none"), srtCol = 45, 
#          offsetCol = 0.1, cexCol = 0.9,
#          ColSideColors = c(
#            rep("red",4),
#            rep("green",4),
#            rep("blue",4))
#)

as.matrix(auc_mtx_mac)

