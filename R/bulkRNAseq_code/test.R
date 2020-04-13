rm(list=ls())
library(readxl)
library("DESeq2")
library("taRifx")
library("BiocParallel")
register(MulticoreParam(4))
library("dplyr")
library("RColorBrewer")
library("gplots")
library("factoextra")
library("vsn")
library("ggplot2")
library("pheatmap")
library("plotly")
library("dichromat")

data <- read.delim("mtx.txt",sep="\t",header=T)
coldata=read.table("coldata.tsv",sep="\t",header=T)
colnames(coldata) <- c("samplename","origin","celltype","subset")
coldata <- coldata[order(coldata$origin, decreasing = TRUE),]
coldata <- remove.factors(coldata)
coldata$celltype <- c(10:41)
coldata <- coldata[order(coldata$subset, coldata$origin, decreasing = FALSE),]
coldata$subset <- factor(coldata$subset)
coldata$origin <- factor(coldata$origin, levels = c("SPF","GF"))
rownames(data) <- data$Transcript

datas <- data[,7:38]
temp <- coldata[order(coldata$samplename, decreasing = FALSE),]
colnames(datas) <- paste0(temp$subset,"_",temp$celltype,"_",temp$origin)
colorder <- colnames(datas)
datas <- datas[,order(colorder, decreasing = FALSE)]
colnames(datas)
nonZeroRows = which(rowSums(datas) > 10)
countData.rounded.nonZero<-datas[names(nonZeroRows), ]

print(paste("There are a total of ",nrow(countData.rounded.nonZero),
            " genes with atleast one sample containing non-zero values out of ",
            nrow(data)," genes"))

#Remove rows with very low variance across all samples (unlikely to be differentially expressed)
getVar <- apply(countData.rounded.nonZero[,-1], 1, var)
quantile(getVar)
param <- 0.5
countData.rounded.nonZero.nonLowVar <- countData.rounded.nonZero[getVar > param & !is.na(getVar), ]

print(paste("There are a total of ",
            nrow(countData.rounded.nonZero.nonLowVar),
            "genes remaining after eliminating genes with low variance"))

boxplot(log2(countData.rounded.nonZero.nonLowVar+1), las=2,
        names=colnames(countData.rounded.nonZero), xlab = "",
        col = c(rep("darkgrey",4),rep("grey",4),rep("blue",4),rep("yellow",4),
                rep("red",4),rep("orange",4),rep("purple",4),rep("green",4)),
        ylab = "Raw read counts per gene (log2)")

countData.rounded.nonZero.nonLowVar=round(countData.rounded.nonZero.nonLowVar,0)
test <- countData.rounded.nonZero.nonLowVar

test <- test[apply(test,1,function(x){(sum(x == 0)/length(x)) < 0.90}),]
test <- test[rowSums(test) >= 10,]

boxplot(log2(test+1), las=2,
        names=colnames(test), xlab = "",
        col = c(rep("darkgrey",4),rep("grey",4),rep("blue",4),rep("yellow",4),
                rep("red",4),rep("orange",4),rep("purple",4),rep("green",4)),
        ylab = "Raw read counts per gene (log2)")

#Construct the data object from the matrix of counts and the metadata table,
ddsMat <- DESeqDataSetFromMatrix(countData = test, 
                                 colData = coldata, design=~subset)
keep <- rowSums(counts(ddsMat)) >= 10
ddsMat <- ddsMat[keep,]
#ddsMat$origin <- relevel(ddsMat$condition, "Naive")
as.data.frame( colData(ddsMat) )
dds <- DESeq(ddsMat)
plotDispEsts(dds)
resultsNames(dds)
res <- results(dds)
plotMA(res, alpha = 0.05, ylim = c(-30,30))

####Index Generation
index <- data[,c("Transcript","Gene.Symbol")]
index <- index[intersect(index$Transcript, rownames(assay(vsd))),]
index <- remove.factors(index)

####DEG Extract
P1toP2 <- results(dds, contrast = c("subset","P1","P2"))
DEG_P1toP2 <- subset(P1toP2, log2FoldChange > 10 & padj < 0.05)
DEG_P1toP2 <- DEG_P1toP2[order(DEG_P1toP2$log2FoldChange, decreasing = TRUE),]
write.table(P1toP2, file = "P1toP2_posonly.txt", sep = '\t')

P2toP3 <- results(dds, contrast = c("subset","P2","P3"))
DEG_P2toP3 <- subset(P2toP3, log2FoldChange > 10 & padj < 0.05)
DEG_P2toP3 <- DEG_P2toP3[order(DEG_P2toP3$log2FoldChange, decreasing = TRUE),]
write.table(P2toP3, file = "P2toP3_posonly.txt", sep = '\t')

P3toP4 <- results(dds, contrast = c("subset","P3","P4"))
DEG_P3toP4 <- subset(P3toP4, log2FoldChange > 10 & padj < 0.05)
DEG_P3toP4 <- DEG_P3toP4[order(DEG_P3toP4$log2FoldChange, decreasing = TRUE),]
write.table(P3toP4, file = "P3toP4_posonly.txt", sep = '\t')

P4toP5 <- results(dds, contrast = c("subset","P4","P5"))
DEG_P4toP5 <- subset(P4toP5, log2FoldChange > 10 & padj < 0.05)
DEG_P4toP5 <- DEG_P4toP5[order(DEG_P4toP5$log2FoldChange, decreasing = TRUE),]
write.table(P4toP5, file = "P4toP5_posonly.txt", sep = '\t')

P5toP6 <- results(dds, contrast = c("subset","P5","P6"))
DEG_P5toP6 <- subset(P5toP6, abs(log2FoldChange) > 10 & padj < 0.05)
DEG_P5toP6 <- DEG_P5toP6[order(DEG_P5toP6$log2FoldChange, decreasing = TRUE),]
write.table(P5toP6, file = "P5toP6_posonly.txt", sep = '\t')

P4toP7 <- results(dds, contrast = c("subset","P4","P7"))
DEG_P4toP7 <- subset(P4toP7, log2FoldChange > 10 & padj < 0.05)
DEG_P4toP7 <- DEG_P4toP7[order(DEG_P4toP7$log2FoldChange, decreasing = TRUE),]
write.table(P4toP7, file = "P4toP7_posonly.txt", sep = '\t')

P7toP8 <- results(dds, contrast = c("subset","P7","P8"))
DEG_P7toP8 <- subset(P7toP8, abs(log2FoldChange) > 10 & padj < 0.05)
DEG_P7toP8 <- DEG_P7toP8[order(DEG_P7toP8$log2FoldChange, decreasing = TRUE),]
write.table(P7toP8, file = "P7toP8_posonly.txt", sep = '\t')

DEGList <-
c(head(rownames(DEG_P1toP2),10), head(rownames(DEG_P2toP3),10), head(rownames(DEG_P3toP4),10),
  head(rownames(DEG_P4toP5),10), head(rownames(DEG_P5toP6),10), head(rownames(DEG_P4toP7),10),
  head(rownames(DEG_P7toP8),10))
DEGList <- DEGList[!duplicated(DEGList)]

plotMA(P1toP2, alpha=0.05, ylim=c(-30,30))

metadata(res)$alpha
metadata(res)$filterThreshold
plot(metadata(P1toP2)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(P1toP2)$lo.fit, col="red")
abline(v=metadata(P1toP2)$filterTheta)

head(assays(dds)[["mu"]])
head(assays(dds)[["cooks"]])
head(dispersions(dds))
head(mcols(dds)$dispersion)

#dds$group <- factor(paste0(dds$origin, dds$subset))
#design(dds) <- ~ group
#dds <- DESeq(dds)
#resultsNames(dds)

qs <- c(0, quantile( P1toP2$baseMean[P1toP2$baseMean > 0], 0:7/7 ))
bins <- cut(P1toP2$baseMean, qs)
levels(bins) <- paste0("~", round(.5*qs[-1]+.5*qs[-length(qs)]))
ratios <- tapply(P1toP2$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE))
barplot(ratios, xlab = "mean normalized count", ylab = "ratio of small $p$ values")

vsd <- varianceStabilizingTransformation(dds)
head(assay(vsd),3)

####SD Plot
msd <- meanSdPlot(assay(vsd))

msd$gg + ggtitle("SdPlot") + 
  scale_fill_gradient(low = "yellow", high = "darkred") + 
  scale_y_continuous(limits = c(0, 8))  

boxplot(log2(assay(vsd)), las=2, names=colnames(countData.rounded.nonZero),
        col = c(rep("darkgrey",4),rep("grey",4),rep("blue",4),rep("yellow",4),
                rep("red",4),rep("orange",4),rep("purple",4),rep("green",4)))

#### Heatmap of the count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("subset","origin")])
gene.sym <- rownames(assay(vsd)[select,])
gene.sym <- index[gene.sym,]$Gene.Symbol
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, labels_row = gene.sym)

#### PCA Plot
Z=plotPCA(vsd, intgroup = c("subset", "origin"), returnData=TRUE)
percentVar <- round(100 * attr(Z, "percentVar"))
ggplot(Z, aes(PC1, PC2, color=subset, shape=origin)) +
  geom_point(size=4) + scale_color_hue(l=55,c=80) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position = "right", legend.box = "horizontal",
        legend.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.text = element_text(size = 20))+
  coord_fixed()
  #geom_label(aes(label = colnames(test)))

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
fviz_eig(pca)
condition=as.vector(as.matrix(coldata$subset))
scores <- data.frame(cbind(pca$x[,1:3],condition))
percentVariance <- round(100 * percentVar[1:3])

## 3D PCA Plot
plot_ly(scores,x=~PC1,y=~PC2,z=~PC3,color=~condition,text=~as.matrix(row.names(scores)))%>%
  add_markers()%>%
  layout(scene = list(xaxis = list(title = paste('PC1:',percentVariance[1],"% variance")),
                      yaxis = list(title = paste('PC2:',percentVariance[2],"% variance")),
                      zaxis = list(title = paste('PC3:',percentVariance[3],"% variance"))))

#### VolcanoPlot
library(tidyr)
library(ggrepel)
library(dplyr)

res <- P1toP2
res$ensembl <- rownames(res)
mtx <- as.data.frame(res) %>% drop_na()
mutateddf <- mutate(mtx, sig=ifelse(-log10(mtx$pvalue)>2 & abs(mtx$log2FoldChange)> 2 & mtx$padj<0.01, "Sig", "N.S.")) #Will have different colors depending on significance
input <- mutateddf#convert the rownames to a column
volc = ggplot(input, aes(log2FoldChange, -log10(padj))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance
  geom_vline(xintercept = c(-2,2), colour = "red", linetype = "dashed")+
  geom_hline(yintercept = 2, colour = "blue", linetype = "dashed")+
  scale_color_manual(values=c("black", "red")) + 
  ggtitle("s1 vs s5") #e.g. 'Volcanoplot DESeq2'
volc
a <- subset(input, input$sig=="Sig"& abs(input$log2FoldChange) > 15)
a$symbol <- index[a$ensembl,]$Gene.Symbol
volc+geom_text_repel(data=a, aes(label=symbol))#adding text for the top 20 genes

## Extracting PC genes
PC1genes <- data.frame(PC1 = sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:51])
PC1genes$symbol <- index[rownames(PC1genes),]$Gene.Symbol

PC2genes <- data.frame(PC2 = sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:63])
PC2genes$symbol <- index[rownames(PC2genes),]$Gene.Symbol

##### PC genes & MP maturation genes & Heatmap
mp.mature <- c("Cd9","Cx3cr1","Cd63","Adgre1","Cd68","Ccr2")
b <- unlist(sapply(mp.mature, function(x) rownames(index[index$Gene.Symbol == x,]) ))
b <- as.character(b)

pc1gene <- c(rownames(PC1genes))
pc2gene <- c(rownames(PC2genes))
pcgene <- setdiff(pc1gene,pc2gene)

# vsd matrix
test <- assay(vsd)
a <- test[b,] # MP maturatoin genes
a <- test[pcgene,] # PC1, 2 genes

#a$gene.symbol <- index[b,]$Gene.Symbol
nonZeroRows = which(rowSums(a) != 0)
a<-a[names(nonZeroRows), ]
a <- as.matrix(a)
rownames(a) <- index[rownames(a),]$Gene.Symbol
#integrate the gene expression value of isoforms
a <- integration(a)

#cData <- t(apply(a,1,function(x)x-median(x)))
scale.range <- c(-4,4)
scale.brakes <- seq(scale.range[1],scale.range[2], by = 0.1)
n.brakes <- length(scale.brakes) - 1
#Palette <- colorpanel(n.brakes,'lightcyan1','lightyellow','indianred2')
Palette <- colorpanel(n.brakes,'darkblue','lightgoldenrodyellow','darkred')
#Palette <- colorpanel(n.brakes,'cyan','red')
require(scales)
my_color_palette <- hue_pal()(8)

rowDendro <- as.dendrogram(hclust(as.dist(1-cor(t(a)))))
#colDendro <- as.dendrogram(hclust(as.dist(1-cor(a))))
heatmap.2(a, trace = 'none', scale = 'row',
          offsetRow = 0.01, 
          cexRow = 1.0,col = Palette, #breaks = scale.brakes, 
          symkey = T,
          Rowv = rowDendro, Colv = 'NA', 
          dendrogram = c("row"), srtCol = 45, 
          offsetCol = 0.1, cexCol = 0.0001,
          ColSideColors = c(rep(my_color_palette[1],4),
                            rep(my_color_palette[2],4),
                            rep(my_color_palette[3],4),
                            rep(my_color_palette[4],4),
                            rep(my_color_palette[5],4),
                            rep(my_color_palette[6],4),
                            rep(my_color_palette[7],4),
                            rep(my_color_palette[8],4))
          )

####Integration Functoin Definition
integration <- function(a){
  newmtx = data.frame()
  temp = data.frame()
  for (i in 1:length(rownames(a))){
    dup <- sum(rownames(a)[i] == rownames(a))
    deci <- (i > grep(rownames(a)[i], rownames(a))[1])
    if (!deci){
     if (dup > 1){
      temp = colSums(a[rownames(a)[i] == rownames(a),])
      temp <- t(as.data.frame(temp))
      rownames(temp) = rownames(a)[i]
      } else {
      temp = a[rownames(a)[i],]
      temp <- t(as.data.frame(temp))
      rownames(temp) = rownames(a)[i]
      }
      newmtx <- rbind(newmtx, temp)
    }
  }
  return(as.matrix(newmtx))
}
