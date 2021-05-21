# READ alignment files ###################################################
# Test the GNN for occurance of PFAMs associated with electron transfer
# The script gives an output of clusters which contain one of these PFAM
# The cutoff for co-occurance in this analysis was 20% meaning only those clusters
# in Which the PFAM occurs in the neighborhood 
library(dplyr)
library(stringr)
library(plyr)
library(ggplot2)
library(ggsci)

dir_name <-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_name)

# Read the files needed
BATSGNN<-read.csv(file = "BATS_1_SSN_150 Full Network GNN default node.csv")
# list of PFAM ids of interest
PFAMs_full<-read.csv(file = "PFAM_ids_ECP.csv", header=TRUE)
PFAMS<-PFAMs_full$PFAM

# Subset the GNN of cluster containing at least one of the PFAMs of interest
BATSGNN_ECP<-BATSGNN[grep(paste(PFAMS, collapse="|"), BATSGNN$Pfam), ]
SNN_numbers<-as.data.frame(table(BATSGNN_FAS$SSN.Cluster.Number))

# We need to separate the PFAMs, because there are rows with multiple in one cell
PFAMs_red<-c()
for (n in 1:nrow(BATSGNN_ECP)){
  newPFAM<-BATSGNN_ECP$Pfam[n]
  if (grepl( "-",BATSGNN_ECP$Pfam[n], fixed = TRUE)){
    newPFAM<-unlist(str_split(BATSGNN_ECP$Pfam[n], "-"))
  }
  PFAMs_red<-c(PFAMs_red,newPFAM)
}
PFAM_counts<-count(data.frame(random=1:length(PFAMs_red),PFAM_ids=PFAMs_red),"PFAM_ids")
PFAM_counts <-PFAM_counts[order(PFAM_counts$freq),]

PFAM_counts <- merge(PFAM_counts, PFAMs_full, by.x = "PFAM_ids", by.y = "PFAM")

p<-ggplot(data=PFAM_counts, aes(x=reorder(PFAM_ids, -freq), y=freq)) +
  geom_bar(stat="identity")+
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(y=20,label=Desc),angle = 90, vjust=0.2, color="black", size=3)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        text = element_text(size = 9, family = "Arial", color = "grey"),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=0.53))+
  scale_fill_npg()
p


length(unique(BATSGNN$Pfam))
length(unique(BATSGNN$SSN.Cluster.Number))

length(unique(BATSGNN_ECP$SSN.Cluster.Number))
sum(PFAM_counts$freq)

#Issue PFAM
# Radical SAM enzymes are often found together
# Radical SAM
PF04055
#ThiH and HydG
PF00037
PF02906
PF10588
PF12838
PF13510
PF02256
PF13353
PF13394

# Fusobacterium cluster 157, 4Fe-4S ferredoxin
# Dialister cluster 146, mulitple biotin synthases
# Desulfoluna cluster 896, same ferredoxin as fusobacterium
# Methylobacterium sp cluster 55, small ferredoxin, ferredoxin reductase
# Chroococcidiopsis. 50
# Desulfurella multipotens, Orenia, small ferredoxin 605, 682
# Veillonellaceae bacterium 146, random sam radical
# BioB with fumurate reductase/SDH 226, 241, 286, 513, 529
# BioB with NADPH dehydrogenases, 354, 800, 903, 1079
# BioB with Epoxyqueuosine reductase, 11, 605, 904



# Clusters with ECPs
ECP_clusters<-c(157, 896, 55, 50, 605, 682)
