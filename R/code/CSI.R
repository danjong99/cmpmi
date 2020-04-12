library(magrittr)
library(tibble)
library(expss)

## read auc matrix
auc_mtx <- read.table("auc_mtx.txt", header=T, row.names=1)
regulons <- rownames(auc_mtx)

##############################################################################

## get cell barcode information
cell_1 <- colnames(auc_mtx)
cellinfo <- read.table("cellInfo.txt", header=T, row.names=1)

head(cellinfo)
cellinfo %>% dplyr::group_by(CellType) %>% dplyr::count()

sorted_colnames <- cellinfo %>% dplyr::arrange(CellType) %>% dplyr::select(CellID) %>% as.matrix() %>% as.vector()
names(sorted_colnames) <- cellinfo %>% dplyr::arrange(CellType) %>% dplyr::select(CellType) %>% as.matrix() %>% as.vector()

temp <- sorted_colnames[sorted_colnames %in% cell_1]
names(temp)

auc_mtx_sorted <- auc_mtx[,temp]
auc_mtx_sorted_renames <- auc_mtx_sorted
colnames(auc_mtx_sorted_renames) <- names(temp)

###################################################################################

rownames(cov(auc_mtx_sorted))
rownames(cov(t(auc_mtx_sorted)))

pcc <- cor(t(auc_mtx_sorted))
pcc[1:5,1:5]

i=0
## return # of element which is lower than ref
count_lower <- function(ref, vect) {
  i=0
  for(elem in vect){if(elem < ref){i=i+1}}
  return(i)
}

csi <- matrix(nrow=215,ncol=215)
rownames(csi) <- colnames(pcc)
colnames(csi) <- colnames(pcc)

for(i in 1:nrow(pcc)){
  for(j in 1:nrow(pcc)){
    if(i==j) {
      csi[i,j]=1
      }
    else {
      csi[i,j]=( (count_lower(pcc[i,j],pcc[i,]) - 0.05) / 215 )
      }
  }
}


library(corrplot)
corrplot(csi, method="color", order = "hclust", tl.cex = 0.1,
         tl.col = "white", cl.lim = c(-0.05,1))
