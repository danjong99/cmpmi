library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
library(GOSemSim)
library(DOSE)
library(RDAVIDWebService)

cmp <- readRDS("seurat_combnied.rds")
C2vsAll <- FindMarkers(object = cmp, ident.1 = 2)
C2vsAll_go <- subset(C2vsAll, abs(avg_logFC) > 0.5)
id1 <- rownames(subset(C2vsAll_go, avg_logFC > 0))
id2 <- rownames(subset(C2vsAll_go, avg_logFC < 0))
C2vsC6 <- FindMarkers(object = cmp, ident.1 = 2, ident.2 = 6) 
C2vsC6_go <- subset(C2vsC6, abs(avg_logFC) > 0.5)
id3 <- rownames(subset(C2vsC6_go, avg_logFC > 0))
id4 <- rownames(subset(C2vsC6_go, avg_logFC < 0))

id1 <- bitr(id1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
id2 <- bitr(id2, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
id3 <- bitr(id3, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
id4 <- bitr(id4, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")

C2vsAll_pos <- enrichGO(gene = id1$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'ENTREZID',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           readable = TRUE)
C2vsAll_neg <- enrichGO(gene = id2$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = 'ENTREZID',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable = TRUE)
C2vsC6_pos <- enrichGO(gene = id3$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'ENTREZID',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           readable = TRUE)
C2vsC6_neg <- enrichGO(gene         = id4$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = 'ENTREZID',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable = TRUE)

C2vsAll_pos <- simplify(C2vsC4_pos, cutoff=0.7, by="p.adjust", select_fun=min)
C2vsAll_neg <- simplify(C2vsC4_neg, cutoff=0.7, by="p.adjust", select_fun=min)
C2vsC6_pos <- simplify(C2vsC6_pos, cutoff=0.7, by="p.adjust", select_fun=min)
C2vsC6_neg <- simplify(C2vsC6_neg, cutoff=0.7, by="p.adjust", select_fun=min)

dotplot(dropGO(C2vsAll_pos, level = 1))
dotplot(dropGO(C2vsAll_neg, level = 1))
dotplot(dropGO(C2vsC6_pos, level = 1))
dotplot(dropGO(C2vsC6_neg, level = 1))

dotplot(C2vsAll_pos)
dotplot(C2vsAll_neg)
dotplot(C2vsC6_pos)
dotplot(C2vsC6_neg)

emapplot(C2vsAll_pos)
emapplot(C2vsAll_neg)
emapplot(C2vsC6_pos)
emapplot(C2vsC6_neg)

write.table(C2vsAll_pos, file = "C2vsAll_Positive.txt", sep = '\t')
write.table(C2vsAll_neg, file = "C2vsAll_Negative.txt", sep = '\t')
write.table(C2vsC6_pos, file = "C2vsC6_Positive.txt", sep = '\t')
write.table(C2vsC6_neg, file = "C2vsC6_Negative.txt", sep = '\t')