###########################################################
### Branched Heatmap & TF Target Genes 
###########################################################

atf3Target0.9 <- atf3TargetsInfo %>% 
  filter(highConfAnnot == TRUE, Genie3Weight > quantile(atf3TargetsInfo$Genie3Weight, 0.9))
a1 <- branchedHeatmap[intersect(c(atf3Target0.9$gene, "Atf3"), rownames(branchedHeatmap)),]

atf4Target0.9 <- atf4TargetsInfo %>% 
  filter(highConfAnnot == TRUE, Genie3Weight > quantile(atf4TargetsInfo$Genie3Weight, 0.9))
a2 <- branchedHeatmap[intersect(c(atf4Target0.9$gene, "Atf4"), rownames(branchedHeatmap)),]

bhlhe41Target0.9 <- Bhlhe41TargetsInfo %>% 
  filter(highConfAnnot == TRUE, Genie3Weight > quantile(Bhlhe41TargetsInfo$Genie3Weight, 0.9))
a3 <- branchedHeatmap[intersect(c(bhlhe41Target0.9$gene, "Bhlhe41"), rownames(branchedHeatmap)),]

bach1Target0.9 <- bach1TargetsInfo %>% filter(Genie3Weight > 0.01)
branchedHeatmap[intersect(c(bach1Target0.9$gene, "Bach1"), rownames(branchedHeatmap)),]

cebpbTarget0.9 <- cebpbTargetsInfo %>% 
  filter(highConfAnnot == TRUE, Genie3Weight > 0.01)
branchedHeatmap[intersect(c(cebpbTarget0.9$gene, "Cebpb"), rownames(branchedHeatmap)),]

cremTarget0.9 <- CremTargetsInfo %>% 
  filter(highConfAnnot == TRUE, Genie3Weight > 0.01)
a6 <- branchedHeatmap[intersect(c(cremTarget0.9$gene, "Crem"), rownames(branchedHeatmap)),]

fosbTarget0.9 <- FosbTargetsInfo %>% 
  filter(highConfAnnot == TRUE, Genie3Weight > 0.01)
branchedHeatmap[intersect(c(fosbTarget0.9$gene, "Fosb"), rownames(branchedHeatmap)),]

fosTarget0.9 <- FosTargetsInfo %>% 
  filter(Genie3Weight > 0.01)
branchedHeatmap[intersect(c(fosTarget0.9$gene, "Fosb"), rownames(branchedHeatmap)),]

junbTarget0.9 <- junbTargetsInfo %>% 
  filter(Genie3Weight > 0.01)
branchedHeatmap[intersect(c(junbTarget0.9$gene, "Junb"), rownames(branchedHeatmap)),]

jundTarget0.9 <- jundTargetsInfo %>% 
  filter(highConfAnnot == TRUE, Genie3Weight > 0.01)
branchedHeatmap[intersect(c(jundTarget0.9$gene, "Jund"), rownames(branchedHeatmap)),]

klf9Target0.9 <- klf9TargetsInfo %>% filter(Genie3Weight > 0.001)
branchedHeatmap[intersect(c(klf9Target0.9$gene, "Klf9"), rownames(branchedHeatmap)),]

mafkTarget0.9 <- mafkTargetsInfo %>% filter(Genie3Weight > 0.001)
branchedHeatmap[intersect(c(mafkTarget0.9$gene, "Mafk"), rownames(branchedHeatmap)),]

mitfTarget0.9 <- MitfTargetsInfo %>% filter(Genie3Weight > 0.003)
branchedHeatmap[intersect(c(mitfTarget0.9$gene, "Mitf"), rownames(branchedHeatmap)),]

nficTarget0.9 <- NficTargetsInfo %>% filter(Genie3Weight > 0.003)
branchedHeatmap[intersect(c(nficTarget0.9$gene, "Nfic"), rownames(branchedHeatmap)),]

nr1h3Target0.9 <- Nr1h3TargetsInfo %>% filter(Genie3Weight > 0.001)
branchedHeatmap[intersect(c(nr1h3Target0.9$gene, "Nr1h3"), rownames(branchedHeatmap)),]

prdm1Target0.9 <- Prdm1TargetsInfo %>% 
  filter(Genie3Weight > 0.001)
branchedHeatmap[intersect(c(prdm1Target0.9$gene, "Prdm1"), rownames(branchedHeatmap)),]

spicTarget0.9 <- SpicTargetsInfo %>% 
  filter(Genie3Weight > 0.001)
branchedHeatmap[intersect(c(spicTarget0.9$gene, "Spic"), rownames(branchedHeatmap)),]

tcf4Target0.9 <- Tcf4TargetsInfo %>% 
  filter(Genie3Weight > 0.005)
branchedHeatmap[intersect(c(tcf4Target0.9$gene, "Tcf4"), rownames(branchedHeatmap)),]

