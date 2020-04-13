#### This code is only running properly in Monocle2
#### This will not be working in Monocle3

library(monocle)
library(Seurat)

############ Monocle2 #####################
cmp_mac <- subset(cmp, idents = c("1","2","3","4","6","7","11"))

#### For trajectory of colon macrophages, only cells with more than 500 nGenes
#### were taken into account for downstream analysis
Monocle_Obj <- SubsetData(object = cmp_mac,
                          cells.use = colnames(cmp_mac@data[,cmp_mac@meta.data$nGene > 500]))
cmac <- Monocle_Obj@raw.data[rownames(Monocle_Obj@data),colnames(Monocle_Obj@data)]

pd_cmac <- data.frame(cellID=colnames(cmac), cell.ident = Monocle_Obj@ident[colnames(cmac)], orig.ident = Monocle_Obj@meta.data$orig.ident)
identical(rownames(pd_cmac), as.character(pd_cmac$cellID))

fd_cmac <- data.frame(gene_short_name=rownames(cmac), totalReads=rowSums(as.matrix(cmac)))
identical(rownames(fd_cmac), as.character(fd_cmac$gene_short_name))

## construct new monocle CellDataSet
pd <- new("AnnotatedDataFrame", data = pd_cmac)
fd <- new("AnnotatedDataFrame", data = fd_cmac)
cMP.monocle <- newCellDataSet(as.matrix(colonMac), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())

cMP.monocle <- estimateSizeFactors(cMP.monocle)
cMP.monocle <- estimateDispersions(cMP.monocle)

pData(cMP.monocle)$Total_mRNAs <- Matrix::colSums(exprs(cMP.monocle))
upper_bound <- 10^(mean(log10(pData(cMP.monocle)$Total_mRNAs)) + 2*sd(log10(pData(cMP.monocle)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cMP.monocle)$Total_mRNAs)) - 2*sd(log10(pData(cMP.monocle)$Total_mRNAs)))

cMP.monocle <- detectGenes(cMP.monocle, min_expr = 0.1)
cMP.monocle <- cMP.monocle[,pData(cMP.monocle)$Total_mRNAs > lower_bound & pData(cMP.monocle)$Total_mRNAs < upper_bound]
fData(cMP.monocle)$use_for_ordering <- fData(cMP.monocle)$num_cells_expressed > 0.05 * ncol(cMP.monocle)
expressed_genes <- row.names(subset(fData(cMP.monocle), num_cells_expressed >= 10))

### Determine Ordering Genes
disp_table <- dispersionTable(cMP.monocle)
ordering_genes <- subset(disp_table,mean_expression > 0.11 & dispersion_empirical > 1*dispersion_fit)$gene_id
length(ordering_genes)

################# Common Functions for Trajectory
cMP.monocle <- setOrderingFilter(cMP.monocle, ordering_genes)
plot_ordering_genes(cMP.monocle)
cMP.monocle <- reduceDimension(cMP.monocle, max_components=2, method = 'DDRTree', norm_method = 'log')
cMP.monocle <- orderCells(cMP.monocle, reverse = TRUE)
cMP.monocle <- orderCells(cMP.monocle, root_state=GM_state(cMP.monocle))

cluster_cols
colVars$CellType[c("1","2","3","4","6","7","11")]


plot_cell_trajectory(cMP.monocle, color_by="State", cell_size=1)+theme(legend.position = "right")
plot_cell_trajectory(cMP.monocle)+theme(legend.position = "right")

png(filename = "Pseudotime_pseudotime.png",
    width = 3000, height = 3000, units = "px", pointsize = 10, res = 600)
plot_cell_trajectory(cMP.monocle, color_by="Pseudotime", cell_size=1, cell_link_size = 2,theta = 90)+
theme(legend.position = c(0.7,0.75), legend.title=element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "ivory",
                                       size = 0.5, linetype = "solid",
                                       colour = "darkblue"))+
  stat_density2d(color='black', h = 4.1, alpha=I(0.3), size=I(0.3))+
  scale_colour_continuous(trans = 'reverse')
dev.off()

png(filename = "Pseudotime_celltype.png",
    width = 3500, height = 3000, units = "px", pointsize = 10, res = 600)
plot_cell_trajectory(cMP.monocle, color_by="cell.ident", cell_size=1, cell_link_size = 2,theta = 90)+
  theme(legend.position = "right", legend.title=element_blank(),
        legend.text = element_text(face = "bold", size = 20),
        legend.background = element_rect(fill = "ivory",
                                         size = 0.5, linetype = "solid",
                                         colour = "darkblue"),
        legend.spacing.y = unit(0.1, 'cm'))+
  scale_color_manual(values = colVars$CellType[c("1","2","3","4","6","7","11")], name = "cell.ident")+
  stat_density2d(color='black', h = 4.1, alpha=I(0.3), size=I(0.3))+
  guides(colour = guide_legend(override.aes = list(size=3),
                               keyheight = 0.6,
                               default.unit = "inch",
                               reverse = TRUE))+facet_wrap(~cell.ident, nrow = 1)
dev.off()

plot_cell_trajectory(cMP.monocle, color_by="State", cell_size=1)+facet_wrap(orig.ident~cell.ident, nrow = 2)
plot_cell_trajectory(cMP.monocle, color_by="State", cell_size=1)+facet_wrap(~orig.ident, nrow = 1)+theme(legend.position = "right")
plot_cell_trajectory(cMP.monocle, color_by = "cell.ident", cell_size = 1, show_tree = FALSE,
                     use_color_gradient = TRUE, markers = "Junb")+facet_wrap(orig.ident~cell.ident, nrow = 2)+theme(legend.position = "right")

## re-order the tree by Pseudotime
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State,pData(cds)$Pseudotime)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else { return (1) }
}

## check gene expression around a branch (or time point)
## BEAM was excuted in server
BEAM_res1 <- BEAM(cMP.monocle, branch_point=1, cores = 16)
BEAM_res2 <- BEAM(cMP.monocle, branch_point=2, cores = 16)
BEAM_res3 <- BEAM(cMP.monocle, branch_point=3, cores = 16)

BEAM_res1 <- BEAM_res1[order(BEAM_res1$qval),]
BEAM_res1 <- BEAM_res1[,c("gene_short_name", "pval", "qval")]

BEAM_res2 <- BEAM_res2[order(BEAM_res2$qval),]
BEAM_res2 <- BEAM_res2[,c("gene_short_name", "pval", "qval")]

BEAM_res3 <- BEAM_res3[order(BEAM_res3$qval),]
BEAM_res3 <- BEAM_res3[,c("gene_short_name", "pval", "qval")]

expressed_genes <- row.names(subset(fData(cMP.monocle), num_cells_expressed >= 10))
cMP.monocle_filtered <- cMP.monocle[expressed_genes,]
mac_genes <- c("Adgre1","Cd63","Ccr2","Cx3cr1")
my_genes <- row.names(subset(fData(cMP.monocle_filtered),gene_short_name %in% mac_genes))
cMP.monocle_subset <- cMP.monocle_filtered[my_genes,]
plot_genes_in_pseudotime(cMP.monocle_subset, color_by = "cell.ident")+ facet_wrap(orig.ident~cell.ident, nrow = 2)

