library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(readxl)

set.seed(345)
TFs_OVERLAP <- read_excel("RSS/TFs_OVERLAP.xlsx", sheet = "RSS.Z")
col_name <- paste0(rep(c("RSSZ","AUC"),13),rep(seq_len(13), each = 2) - 1)
col_name <- paste0(rep(c("RSSZ","AUC"),13))
colnames(TFs_OVERLAP) <- c("regulon",col_name)
TFs_OVERLAP <- as.data.frame(TFs_OVERLAP)
rownames(TFs_OVERLAP) <- TFs_OVERLAP$regulon
TFs_OVERLAP$regulon <- NULL
TFs_OVERLAP1 <- as.data.frame(TFs_OVERLAP)
TFs_OVERLAP2 <- as.data.frame(TFs_OVERLAP)

par(mfrow=c(3,5))
for (i in 1:13) {
  plot(TFs_OVERLAP1[,(i*2-1):(i*2)], type = "p", col = "red",
       lwd = 2, main = paste("Cluster:",c(i-1)), 
       xlim = range(-4:4),
       xlab = "RSSZ", ylab = "AUC",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       pch = 19)
  abline(h=0.2, v=1.0)
}

dev.off()
a <- TFs_OVERLAP1[,1:2]
with(a, symbols(x=a$RSSZ0, y=a$AUC0, circles=sqrt(a$AUC0/pi),inches=1/10,
                  ann=F, bg="red", fg=NULL),
     lwd = 2, main = paste("Cluster:",c(i-1)), 
     xlim = range(-4:4),
     xlab = "RSSZ", ylab = "AUC",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     pch = 19)
abline(h=0.2, v=1.0)


r <- data.frame()
a <- data.frame()
for (i in 1:13) {
a <- TFs_OVERLAP2[,c(1,(i*2):(i*2+1))] %>% filter(RSSZ > 1.0 & AUC > 0.2)
if(nrow(a) != 0){
a$cluster <- c(i-1)
}
r <- rbind(r,a)
}

write.table(r, file = "TFs_OVERLAP.txt", sep = '\t')
a <- TFs_OVERLAP[,c(1,4:5)] %>% filter(RSSZ > 1.0 & AUC > 0.2)
nrow(a)

rss_auc_table <- list()
for (i in 1:13) {
temp <- TFs_OVERLAP[,(i*2-1):(i*2)]
#temp <- TFs_OVERLAP[,1:2]
colnames(temp) <- c("RSSZ","AUC")
temp$regulons <- rownames(temp)
rownames(temp) <- NULL
temp <- as.data.frame(temp)

temp <- mutate(temp, sig=ifelse((temp[,1] > 1.0), "Sig", "N.S."))
temp <- mutate(temp, sig= apply(temp,1,function(x){
                  #x <- unname(x)
                  if(x[4] != "Sig"){
                  ifelse((x[2] > 0.3),"Comm","N.S.")
                  }
                  else return("Sig")
                  }))
rss_auc_table[[i]] <- temp
}

for (i in 1:13) {
  temp = rss_auc_table[[i]]
  aa <- subset(temp, temp$sig==c("Sig")|temp$sig==c("Comm"))
  aa <- rbind(aa %>% filter(AUC > 0.3), aa %>% filter(RSSZ > 1))
  aa <- aa[order(aa$RSSZ, decreasing = TRUE),]
  x_limits <- c(NA, 1.5)
  print(paste("Cluster:",i-1))
  print(aa)
}

for (i in 1:2) {
  temp = rss_auc_table[[7]]
  volc = ggplot(temp, aes(RSSZ, AUC, color = factor(sig))) +
    geom_point(size = 4) +
    #scale_color_manual(values=c("blue","black","red")) + 
    ggtitle(paste("Cluster:",6,"Regulons")) +
    geom_hline(yintercept=0.3, linetype="dashed", color = "red") +
    geom_vline(xintercept=1.0, linetype="dashed", color = "red")
  
  aa <- subset(temp, temp$sig==c("Sig")|temp$sig==c("Comm"))
  aa <- rbind(aa %>% filter(AUC > 0.35), aa %>% filter(RSSZ > 1.5))
  aa <- aa[!duplicated(aa$regulons),]
  aa <- aa[order(aa$RSSZ, decreasing = TRUE),]
  x_limits <- c(NA, 1.5)
  
  p1 <- volc+geom_text_repel(data=aa, size=3, color = "black", aes(label=regulons),
                       arrow = arrow(length = unit(0.03, "npc"), 
                                     type = "closed", ends = "first"),
                       force = 10, xlim = x_limits)+theme_bw()
  grid.arrange(
    tableGrob(aa),
    p1,
    ncol = 2,
    widths = c(1,1),
    clip = FALSE
  )
}



