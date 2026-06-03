rm(list=ls()) 
setwd("/Users/starp/Documents/Dropbox/RNAseq_stuff/topGO")

BiocManager::install("topGO")
BiocManager::install("Rgraphviz")
library(topGO)
library(Rgraphviz)
library(ggplot2)

#This script was originally written by Amanda Katzer (for her CRABS-CLAW paper) and modified by THD for my data! 
#I also used this online tutorial for some info while annotaing this script: https://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

#Call in GO annotations for each of the genes
#I obtained this file from Amanda
#first column contains P.bar gene IDs (ie: Pbar_2022_maker4_00001). Luckily, these gene IDs are present in P.smallii annotation and I have pulled them into my DE files in a previous step with python script on cluster.
#second column indicates the GO terms associated with that gene
geneID2GO <- readMappings(file="./arabidopsis_GO_terms.txt")

#Make a list of the gene names
geneUniverse <- names(geneID2GO)

#########################
## Stage 1 hirs D vs V ##
#########################
#Read in DE data that contains the P.bar gene ID
DEseq_genes <- read.csv("./", header = T)
sub_genes <- DEseq_genes[DEseq_genes$padj < 0.01,]

genesOfInterest <- sub_genes

#I used mRNA_IDs here because that is the column of info that corresponds with the gene IDs in the readMappings file
genesOfInterest <- as.character(genesOfInterest$mRNA_IDs) 
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#Explanation of the "ontology arguments" in topGO from: doi: 10.1007/978-1-4939-3743-1_2
#a gene encodes a gene product, and that gene product carries out a molecular-level process or activity (molecular function) 
#in a specific location relative to the cell (cellular component)
#and this molecular process contributes to a larger biological objective (biological process) comprised of multiple molecular-level processes.

#putting data together into an R object - Biological Processes 
myGOdata <- new("topGOdata", description="stage 1 hirs D vs V, biological processes", 
                ontology="BP", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)

#Fisher's exact test is used to determine if there are nonrandom associations between two categorical variables
#Note: the "weight01" option takes the GO hierarchy into account which is more desirable than the "classic" option which tests each GO term independently 
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO <- usedGO(object=myGOdata)
allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultsFisher", ranksOf="classicFisher", topNodes=length(allGO))
allRes['p_adj'] <- c(allRes$classicFisher, method="fdr")
#list the top ten significant results found
top10Res <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)

#save this data set on local computer
write.csv(allRes, file="./Stg1hirsDvsV_BioP.csv")

#Visualize the position of the statistically significant GO terms in the GO hierarchy
#The significant GO terms are shown as rectangles in the picture. 
#The most significant terms are coloured red and least significant in yellow
par(cex = 0.25) #text size in subsequent plot
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 3, useInfo = 'all')

#putting data together into an R object - Molecular Function
myGOdata <- new("topGOdata", description="stage 1 hirs D vs V, molecular function", 
                ontology="MF", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO <- usedGO(object=myGOdata)
allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultsFisher", ranksOf="classicFisher", topNodes=length(allGO))
allRes['p_adj'] <- p.adjust(allRes$classicFisher, method="fdr")
#list the top ten significant results found
top10Res <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)

#save this data set on local computer
write.csv(allRes, file="./Stg1hirsDvsV_MolF.csv")

par(cex = 0.25) #text size in subsequent plot
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')

#putting data together into an R object - Cellular Component
myGOdata <- new("topGOdata", description="stage 1 hirs D vs V, cellular component", 
                ontology="CC", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO <- usedGO(object=myGOdata)
allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultsFisher", ranksOf="classicFisher", topNodes=length(allGO))
allRes['p_adj'] <- p.adjust(allRes$classicFisher, method="fdr")
#list the top ten significant results found
top10Res <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)

#save this data set on local computer
write.csv(allRes, file="./Stg1hirsDvsV_CellC.csv")

par(cex = 0.25) #text size in subsequent plot
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')

#########################
## Stage 2 hirs D vs V ##
#########################
#Read in DE data 
DEseq_genes <- read.csv("../DESeq2Out/OutputWithAthalInfo/NEW3.6.hirsStg2D_vs_hirsStg2V_with_Athal.csv", header = T)
sub_genes <- DEseq_genes[DEseq_genes$padj < 0.01,]

genesOfInterest <- sub_genes

genesOfInterest <- as.character(genesOfInterest$mRNA_IDs)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

#putting data together into an R object - Biological Processes 
myGOdata <- new("topGOdata", description="stage 2 hirs D vs V, biological processes", 
                ontology="BP", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO <- usedGO(object=myGOdata)
allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultsFisher", ranksOf="classicFisher", topNodes=length(allGO))
allRes['p_adj'] <- p.adjust(allRes$classicFisher, method="fdr")

write.csv(allRes, file="./Stg2hirsDvsV_BioP.csv")

par(cex = 0.25) #text size in subsequent plot
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 16, useInfo = 'all')

#putting data together into an R object - Molecular Function
myGOdata <- new("topGOdata", description="stage 2 hirs D vs V, molecular function", 
                ontology="MF", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO <- usedGO(object=myGOdata)
allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultsFisher", ranksOf="classicFisher", topNodes=length(allGO))
allRes['p_adj'] <- p.adjust(allRes$classicFisher, method="fdr")

write.csv(allRes, file="./Stg2hirsDvsV_MolF.csv")

par(cex = 0.25) #text size in subsequent plot
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'all')

#putting data together into an R object - Cellular Component
myGOdata <- new("topGOdata", description="stage 2 hirs D vs V, cellular component", 
                ontology="CC", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO <- usedGO(object=myGOdata)
allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultsFisher", ranksOf="classicFisher", topNodes=length(allGO))
allRes['p_adj'] <- p.adjust(allRes$classicFisher, method="fdr")

write.csv(allRes, file="./Stg2hirsDvsV_CellC.csv")

par(cex = 0.25) #text size in subsequent plot
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 3, useInfo = 'all')

###########################
# GO-Term Figure creation #
###########################
#read in data frame that contains significantly enriched genes in both early and late P. hirsutus D vs V comparison
TopGoSignif <- read.csv("Top25GOhirsDvsVbothStages.csv") #20568 genes
colnames(TopGoSignif)

TopGo25 <- ggplot(TopGoSignif, aes(x = Stage, y = p.value)) + 
  geom_jitter(aes(color = Significant), size = 4, alpha = 0.85) +
  geom_hline(yintercept = 0.05, color = "red", size = 0.7, linetype = "dashed") +
  xlab("") + ylab("p-value") + ggtitle("Top 25 GO-Terms") +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black"), 
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(size = 0.25, linetype = "dashed", colour = "gray70"), 
        panel.grid.minor = element_line(size = 0.25, linetype = "dashed", colour = "gray70")) +
  scale_x_discrete(limits = c("Early", "Late")) + 
  scale_color_manual(values = c("grey", "blue"))
plot(TopGo25) 


