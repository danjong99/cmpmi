############ Monocle2 #####################
Monocle_Obj <- SubsetData(object = Exr_mac, cells.use = colnames(Exr_mac@data[,Exr_mac@meta.data$nGene > 500]))
colonMac <- Monocle_Obj@raw.data[rownames(Monocle_Obj@data),colnames(Monocle_Obj@data)]

pd_colonMac <- data.frame(cellID=colnames(colonMac), cell.ident = Monocle_Obj@ident[colnames(colonMac)], orig.ident = Monocle_Obj@meta.data$orig.ident)
identical(rownames(pd_colonMac), as.character(pd_colonMac$cellID))

fd_colonMac <- data.frame(gene_short_name=rownames(colonMac), totalReads=rowSums(as.matrix(colonMac)))
identical(rownames(fd_colonMac), as.character(fd_colonMac$gene_short_name))

## construct new monocle CellDataSet
pd <- new("AnnotatedDataFrame", data = pd_colonMac)
fd <- new("AnnotatedDataFrame", data = fd_colonMac)
cMP.monocle <- newCellDataSet(as.matrix(colonMac), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
cMP.monocle <- estimateSizeFactors(cMP.monocle)
cMP.monocle <- estimateDispersions(cMP.monocle)
sizeFactors(cMP.monocle)

pData(cMP.monocle)$Total_mRNAs <- Matrix::colSums(exprs(cMP.monocle))
upper_bound <- 10^(mean(log10(pData(cMP.monocle)$Total_mRNAs)) + 2*sd(log10(pData(cMP.monocle)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cMP.monocle)$Total_mRNAs)) - 2*sd(log10(pData(cMP.monocle)$Total_mRNAs)))

cMP.monocle <- detectGenes(cMP.monocle, min_expr = 0.1)
cMP.monocle <- cMP.monocle[,pData(cMP.monocle)$Total_mRNAs > lower_bound & pData(cMP.monocle)$Total_mRNAs < upper_bound]
fData(cMP.monocle)$use_for_ordering <- fData(cMP.monocle)$num_cells_expressed > 0.05 * ncol(cMP.monocle)
expressed_genes <- row.names(subset(fData(cMP.monocle), num_cells_expressed >= 10))

### Determine Ordering Genes - Unsupervised fashion
disp_table <- dispersionTable(cMP.monocle)
ordering_genes <- subset(disp_table,mean_expression > 0.11 & dispersion_empirical > 1*dispersion_fit)$gene_id
length(ordering_genes)

################# Common Functions for Trajectory
cMP.monocle <- setOrderingFilter(cMP.monocle, ordering_genes)
plot_ordering_genes(cMP.monocle)
cMP.monocle <- reduceDimension(cMP.monocle, max_components=2, method = 'DDRTree', norm_method = 'log')
cMP.monocle <- orderCells(cMP.monocle, reverse = TRUE)
cMP.monocle <- orderCells(cMP.monocle, root_state=GM_state(cMP.monocle))

plot_cell_trajectory(cMP.monocle, color_by="State", cell_size=1)+theme(legend.position = "right")
plot_cell_trajectory(cMP.monocle, color_by="Pseudotime", cell_size=1)+theme(legend.position = "top")
plot_cell_trajectory(cMP.monocle, color_by="cell.ident", cell_size=1)+theme(legend.position = "right")#+facet_wrap(category~.)

plot_cell_trajectory(cMP.monocle, color_by="State", cell_size=1)+facet_wrap(orig.ident~cell.ident, nrow = 2)
plot_cell_trajectory(cMP.monocle, color_by="State", cell_size=1)+facet_wrap(~orig.ident, nrow = 1)+theme(legend.position = "right")
plot_cell_trajectory(cMP.monocle, color_by = "cell.ident", cell_size = 1)+facet_wrap(~cell.ident, nrow = 2)+theme(legend.position = "right")

## re-order the tree by Pseudotime
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State,pData(cds)$Pseudotime)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else { return (1) }
}

## check gene expression around a branch (or time point)
BEAM_res2 <- BEAM(cMP.monocle, branch_point=2, cores = 16)

BEAM_res2 <- BEAM_res2[order(BEAM_res2$qval),]
BEAM_res2 <- BEAM_res2[,c("gene_short_name", "pval", "qval")]

expressed_genes <- row.names(subset(fData(cMP.monocle), num_cells_expressed >= 10))
cMP.monocle_filtered <- cMP.monocle[expressed_genes,]
my_genes <- row.names(subset(fData(cMP.monocle_filtered),gene_short_name %in% c("Ccl6")))
cMP.monocle_subset <- cMP.monocle_filtered[my_genes,]
plot_genes_in_pseudotime(cMP.monocle_subset, color_by = "cell.ident")+ facet_wrap(orig.ident~.)

WT_Cells <- row.names(subset(pData(cMP.monocle_subset),orig.ident %in% c("WT")))
GF_Cells <- row.names(subset(pData(cMP.monocle_subset),orig.ident %in% c("GF")))
WT_subset <- cMP.monocle_subset[,WT_Cells]
plot_genes_in_pseudotime(WT_subset, color_by = "cell.ident")
GF_subset <- cMP.monocle_subset[,GF_Cells]
plot_genes_in_pseudotime(GF_subset, color_by = "cell.ident")


marker_genes <- rownames(subset(fData(cMP.monocle), gene_short_name %in% c("Ccr2","Cx3cr1","Maf","Adgre1",
                                                                           "Cd9","Cd36","Cd63","Cd68","Mrc1",
                                                                           "Cd80","Cd86","Csf2ra","Csf2rb",
                                                                           "Csf1r","Retnla","Cd163","Atf3",
                                                                           "Hes1","Clec7a","Clec4n","F13a1",
                                                                           "Ccl2","Ccl7") ))
diff_test_res <- differentialGeneTest(cMP.monocle[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cMP.monocle[sig_gene_names,],
                        num_clusters = 3,
                        cores = 4,
                        show_rownames = T)



my_genes <- row.names(subset(fData(cMP.monocle),gene_short_name %in% marker_genes[11:20]))
plot_genes_branched_pseudotime(cMP.monocle[my_genes,],
                               branch_point = 2,
                               color_by = "cell.ident",
                               ncol = 1)


cMP.monocle_w8 <- cMP.monocle
save(cMP.monocle_w8, file = "WT2_GF_combined_monocle_mac_w8.rds")
