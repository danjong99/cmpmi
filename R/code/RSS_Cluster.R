library(magrittr)
library(tibble)

## read auc matrix
auc_mtx <- read.table("auc_mtx.txt", header=T, row.names=1)
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

####################################################

auc_mtx_sorted_Mat <- as.matrix(auc_mtx_sorted_renames)
  
## calculate Pr
Pr <- auc_mtx_sorted_Mat / rowSums(auc_mtx_sorted_Mat)
dim(Pr)
## calculate Pc
nCell <- sort(unname(table(colnames(auc_mtx_sorted_Mat))),decreasing = T)
Pc <- matrix(nrow=13, ncol=464)

for(i in 1:length(nCell)){
  if(i==1){
    a=0; b=nCell[i]; c=464-a-b
  }
  else{
    a=sum(nCell[1:(i-1)]); b=nCell[i]; c=464-a-b
  }
  Pc[i,] <- c(rep(0,a), rep(1,b), rep(0,c))
}
Pc <- Pc/rowSums(Pc)
dim(Pc)

#########################################################################################
## Entropy (https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r)
H <- function(v) {
  v <- v[v > 0]
  return(sum(-v * log(v)))
}

JSD_ty <- function(w,m){
  H((w+m)/2) - (H(w)+H(m))/2
}

## RSS
RSS <- function(R, C) {
  return(1-sqrt(JSD_ty(R,C)))
}
#########################################################################################

Pr <- as.matrix(Pr)
dim(Pr)
Pc <- as.matrix(Pc)
dim(Pc)

result <- matrix(nrow=215, ncol=13)
colnames(result) <- 0:12
rownames(result) <- rownames(Pr)

for(i in 1:215){
  for(j in 1:13){
    m <- Pc[j,]
    w <- Pr[i,]
    result[i,j] <- RSS(m,w)
  }
}

write.csv(result, "RSS.txt", quote = F)
#########################################################################################