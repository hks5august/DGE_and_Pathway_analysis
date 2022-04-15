#Load libraries
library(gage)
library(gageData)
library(pathview)

#Load DESeq2 results file as input file
result <- read.table("EB_4D_vs_2D.txt", sep="\t", header = TRUE)

#Check top 10 lines in results 
head(result)

#list of KEGG pathways genes
data("kegg.sets.mm")

#Index of numbers of pathway
data(sigmet.idx.mm)

#Extract  cleaner gene set of important pathways
kegg.sets.mm <- kegg.sets.mm[sigmet.idx.mm]

#Extract the fold changes from the results
foldchanges <- result$logFC

#Check by head
head(foldchanges)

#Extract gene ids
names(foldchanges) <- result$EntrezID

#Check top 10 values
head(foldchanges)

#Get all KEGG pathways associated with genes
kegg_result <- gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)

#Write results of all pathways in a file
write.table(kegg_result, file="all_pathways.txt", sep="\t", row.names=TRUE,quote=FALSE)

#Check attributes
attributes(kegg_result)

#Check upregulated pathways
head(kegg_result$greater)

#check down regulated pathways
head(kegg_result$less)

#Extract all upregulated pathways data
Upreg_pathways <- data.frame(Pathway_id=rownames(kegg_result$greater), kegg_result$greater)

#write all upregulated pathways into a file
write.table(Upreg_pathways, file="upregulated_pathways.txt", sep="\t", row.names=F, quote=FALSE)

#Extract all downregulated pathways data
Downreg_pathways <- data.frame(Pathway_id=rownames(kegg_result$less), kegg_result$less)

#write all downregulated pathways into a file
write.table(Downreg_pathways, file="Downregulated_pathways.txt", sep="\t", row.names=F, quote=FALSE)

#Visualize one upregulated pathway
pathview(gene.data=foldchanges, pathway.id="00900 ", kegg.native=FALSE, species =  "Mus musculus", same.layer = F)

#Obtain a different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="00900 ", kegg.native=TRUE,species =  "Mus musculus", same.layer = F)

#Extract the top 5 upregulated pathways 
keggrespathways_up <- rownames(kegg_result$greater)[1:5]

#Extract the IDs part of each string
keggresids_up <- substr(keggrespathways_up, start=1, stop=8)
head(keggresids_up)

#write Ids in a file
write.table(keggresids_up, file="top5_upreg_pathways_ids.txt", sep="\t", row.names=F, quote=FALSE)

#Draw plots for top 5 upregulated pathways
pathview(gene.data=foldchanges, pathway.id=keggresids_up, species="hsa")

#Extract top 5 down regulated pathways 
keggrespathways_down <- rownames(kegg_result$less)[1:5]

#Extract the IDs 
keggresids_down = substr(keggrespathways_down, start=1, stop=8)
head(keggresids_down)
#write Ids in a file
write.table(keggresids_down, file="top5_downreg_pathways_ids.txt", sep="\t", row.names=F, quote=FALSE)

#Draw plots for top 5 downregulated pathways
pathview(gene.data=foldchanges, pathway.id=keggresids_down, species="hsa")

#Print keggresids_down
print(keggresids_down)