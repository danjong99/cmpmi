regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")

##Regulons and Targets Commonly Strong in colon Macrophages
FosTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Fos",]
FosbTargetsInfo <- FosbTargetsInfo[!is.na(FosbTargetsInfo$Genie3Weight),]
FosbTarget0.9 <- FosbTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

FosbTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Fosb",]
FosbTargetsInfo <- FosbTargetsInfo[!is.na(FosbTargetsInfo$Genie3Weight),]
FosbTarget0.9 <- FosbTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

Tcf4TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Tcf4",]
Tcf4TargetsInfo <- Tcf4TargetsInfo[!is.na(Tcf4TargetsInfo$Genie3Weight),]
Tcf4Target0.9 <- Tcf4TargetsInfo %>% 
  filter(Genie3Weight > 0.01)

NficTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Nfic",]
NficTargetsInfo <- NficTargetsInfo[!is.na(NficTargetsInfo$Genie3Weight),]
NficTarget0.9 <- NficTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

SpicTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Spic",]
SpicTargetsInfo <- SpicTargetsInfo[!is.na(SpicTargetsInfo$Genie3Weight),]
SpicTarget0.9 <- SpicTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

MitfTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Mitf",]
MitfTargetsInfo <- MitfTargetsInfo[!is.na(MitfTargetsInfo$Genie3Weight),]
MitfTarget0.9 <- MitfTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

Bhlhe41TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Bhlhe41",]
Bhlhe41TargetsInfo <- Bhlhe41TargetsInfo[!is.na(Bhlhe41TargetsInfo$Genie3Weight),]
Bhlhe41Target0.9 <- Bhlhe41TargetsInfo %>% 
  filter(highConfAnnot == TRUE)

Nr1h3TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Nr1h3",]
Nr1h3TargetsInfo <- Nr1h3TargetsInfo[!is.na(Nr1h3TargetsInfo$Genie3Weight),]
Nr1h3Target0.9 <- Nr1h3TargetsInfo %>% 
  filter(Genie3Weight > 0.01)

##Regulons and Targets Particularly strong in Cluster 6
atf3TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Atf3",]
atf3TargetsInfo <- atf3TargetsInfo[!is.na(atf3TargetsInfo$Genie3Weight),]
atf3Target0.9 <- atf3TargetsInfo %>% 
  filter(Genie3Weight > 0.01)

bach1TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Bach1",]
bach1TargetsInfo <- bach1TargetsInfo[!is.na(bach1TargetsInfo$Genie3Weight),]
bach1Target0.9 <- bach1TargetsInfo %>% 
  filter(Genie3Weight > 0.01)

mafkTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Mafk",]
mafkTargetsInfo <- mafkTargetsInfo[!is.na(mafkTargetsInfo$Genie3Weight),]
mafkTarget0.9 <- mafkTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

atf4TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Atf4",]
atf4TargetsInfo <- atf4TargetsInfo[!is.na(atf4TargetsInfo$Genie3Weight),]
atf4Target0.9 <- atf4TargetsInfo %>% 
  filter(Genie3Weight > 0.01)

cebpbTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Cebpb",]
cebpbTargetsInfo <- cebpbTargetsInfo[!is.na(cebpbTargetsInfo$Genie3Weight),]
cebpbTarget0.9 <- cebpbTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

jundTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Jund",]
jundTargetsInfo <- jundTargetsInfo[!is.na(jundTargetsInfo$Genie3Weight),]
jundTarget0.9 <- jundTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

klf9TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Klf9",]
klf9TargetsInfo <- klf9TargetsInfo[!is.na(klf9TargetsInfo$Genie3Weight),]
klf9Target0.9 <- klf9TargetsInfo %>% 
  filter(Genie3Weight > 0.01)

junbTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Junb",]
junbTargetsInfo <- junbTargetsInfo[!is.na(junbTargetsInfo$Genie3Weight),]
junbTarget0.9 <- junbTargetsInfo %>% 
  filter(Genie3Weight > 0.01)

##Regulons and Targets Particularly strong in Cluster 4
Prdm1TargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Prdm1",]
Prdm1TargetsInfo <- Prdm1TargetsInfo[!is.na(Prdm1TargetsInfo$Genie3Weight),]
Prdm1Target0.9 <- Prdm1TargetsInfo %>% 
  filter(Genie3Weight > 0.01)

CremTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$TF == "Crem",]
CremTargetsInfo <- CremTargetsInfo[!is.na(CremTargetsInfo$Genie3Weight),]
CremTarget0.9 <- CremTargetsInfo %>% 
  filter(highConfAnnot == TRUE)

