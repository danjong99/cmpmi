library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
library("tidyverse")

cmp <- readRDS("seurat_combnied.rds")
C2vsAll <- FindMarkers(object = cmp, ident.1 = 2)
C2vsC6 <- FindMarkers(object = cmp, ident.1 = 2, ident.2 = 6) 

C2vsAll$gene <- rownames(C2vsAll)
mutateddf <- mutate(C2vsAll, sig=ifelse(-log10(C2vsAll$p_val)>5 & abs(C2vsAll$avg_logFC)> 0.5 & C2vsAll$p_val_adj<0.05, "Sig", "N.S.")) #Will have different colors depending on significance
input <- mutateddf#convert the rownames to a column
volc = ggplot(input, aes(avg_logFC, -log10(p_val_adj))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance
  geom_vline(xintercept = c(-0.5,0.5), colour = "red", linetype = "dashed")+
  geom_hline(yintercept = 1, colour = "blue", linetype = "dashed")+
  scale_color_manual(values=c("black", "red"))+
  scale_x_continuous(breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),
                              limits=c(-1.5,1.5)) +
  ggtitle("Cluster 2 vs All the other MPs") #e.g. 'Volcanoplot DESeq2'
volc+geom_text_repel(data=subset(input, input$sig=="Sig"), aes(label=gene))#adding text for the top 20 genes

C2vsC6$gene <- rownames(C2vsC6)
mutateddf <- mutate(C2vsC6, sig=ifelse(-log10(C2vsC6$p_val)>5 & abs(C2vsC6$avg_logFC)> 0.5 & C2vsC6$p_val_adj<0.05, "Sig", "N.S.")) #Will have different colors depending on significance
input <- mutateddf#convert the rownames to a column
volc = ggplot(input, aes(avg_logFC, -log10(p_val_adj))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance
  geom_vline(xintercept = c(-0.5,0.5), colour = "red", linetype = "dashed")+
  scale_color_manual(values=c("black", "red")) + 
  ggtitle("Cluster 3 vs 6 Genes") #e.g. 'Volcanoplot DESeq2'
volc+geom_text_repel(data=subset(input, input$sig=="Sig"), aes(label=gene))#adding text for the top 20 genes
