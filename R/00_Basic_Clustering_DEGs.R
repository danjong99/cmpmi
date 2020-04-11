library(Seurat)
library(dplyr)

## Analyzed in Seurat 2.0
set.seed(01234)
## Calling All Three Gene Matrix of DropSeq
Exr1.data=read.table(file = paste0("GSE137927_spf_colon_phagocytes_Dropseq_UMI_raw_counts.txt"),sep="\t",header=TRUE,row.names=1)
Exr2.data=read.table(file = paste0("GSE137927_gf_colon_phagocytes_Dropseq_UMI_raw_counts.txt"),sep="\t",header=TRUE,row.names=1)

## Create seurat objects from counting matrix
Exr1 <- CreateSeuratObject(raw.data = Exr1.data, min.cells = 3, min.genes = 200, project = "spf")
Exr2 <- CreateSeuratObject(raw.data = Exr2.data, min.cells = 3, min.genes = 200, project = "gf")

cmp <- MergeSeurat(object1 = Exr1, object2 = Exr2, add.cell.id1 = "spf", 
                            add.cell.id2 = "gf", project = "spfgf")
rm(Exr1.data, Exr2.data, Exr1, Exr2)

cmp <- FilterCells(object = cmp, subset.names = c("nGene"),low.thresholds = c(200), high.thresholds = c(2500))
cmp <- NormalizeData(object = cmp, normalization.method = "LogNormalize", scale.factor = 10000)
cmp <- FindVariableGenes(object = cmp, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
cmp <- ScaleData(object = cmp, vars.to.regress = c("nUMI"))
GenePlot(object = cmp, gene1="nUMI", gene2 = "nGene")
VlnPlot(object = cmp, features.plot = c("nGene", "nUMI"), nCol=2)
cmp <- RunPCA(object = cmp, pc.genes = cmp@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = cmp, pcs.use = 1:3)
PCAPlot(object = cmp)
PCHeatmap(object = cmp_mac, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
cmp <- ProjectPCA(object = cmp, do.print = FALSE)
cmp <- JackStraw(object = cmp, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = cmp, PCs = 1:15)
PCElbowPlot(object = cmp)
cmp <- RunTSNE(object=cmp, dims.use = 1:11, do.fast = TRUE)
cmp <- FindClusters(object = cmp, reduction.type="pca", dims.use = 1:11, resolution = 1.2)
TSNEPlot(object = cmp, reduction.type = "pca", do.return = TRUE, pt.size = 1.5, do.label=TRUE, label.size = 5)
saveRDS(object = cmp, file = "spf_gf_combined_seurat.rds")

#### Added upon Seurat Upgrade to V3.0
cmp <- readRDS("spf_gf_combined_seurat.rds")
cmp <- UpdateSeuratObject(cmp)
#Figure1a
DimPlot(object = cmp, reduction = 'tsne', pt.size = 0.8, label = T, label.size = 6)

#Figure1b
cmp.markers <- readRDS("cmp.all.markers.rds")
top30 <- cmp.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
DoHeatmap(object = cmp, features = top30$gene)

#Figure1c
genestodisplay <- c("Adgre1","Cd63","Zeb2","Cd68","Fcgr1","Cx3cr1", #macrophage markers
                    "Itgax","Cd24a","Zbtb46","Kit","Flt3","Itgae") #dendritic cell markers
FeaturePlot(object = cmp, features = genestodisplay, min.cutoff = 'q05', max.cutoff = 'q95',
            cols = c("lightgrey","red"))

#Figure2b
cmp@meta.data$flora <- cmp@meta.data$orig.ident
index <- cmp@meta.data$orig.ident == "WT"
cmp@meta.data$flora[index] <- c("SPF")
cmp@meta.data$flora <- factor(cmp@meta.data$flora, levels = c("SPF","GF"))
DimPlot(object = cmp, reduction = 'tsne', pt.size = 0.8, label = T,
        label.size = 6, split.by = 'flora')

#Figure2d
cmp_subsets <- subset(cmp, idents = c(2,4,6))
cmp_subsets.markers <- FindAllMarkers(object = cmp_subsets, only.pos = TRUE)
top40 <- cmp_subsets.markers %>% group_by(cluster) %>% top_n(40, avg_logFC)
avg.cmp_subsets <- AverageExpression(cmp_subsets, assays = "RNA", return.seurat = TRUE)
DoHeatmap(object = avg.cmp_subsets, features = top40$gene, draw.lines = FALSE) + NoLegend()

#Figure2e
FeaturePlot(object = cmp, features = c("Mrc1","Il1r2"), split.by = "flora",
            min.cutoff = 'q05', max.cutoff = 'q95',
            cols = c("lightyellow","brown"))

