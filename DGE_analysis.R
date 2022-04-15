#Load packages
library(DESeq2)
library(EnhancedVolcano)

#Load data
PDXdataset = read.table('https://drive.google.com/uc?export=download&id=1tb_qnLwjeAf0B5gWmi7_6SGOqy44dU8x', header = TRUE,row.names = 1, sep='\t')


#Define data and remove id column
data <- (PDXdataset)

colnames(PDXdataset) -> names
meta_classes<- c('ER','ER','ER','ER', 'ER', 'ER', 'ER', 'TN','TN','TN','TN','TN', 'TN', 'TN')
meta_data <- data.frame(names, meta_classes)


#Perform first step of DESeq (Define DESeq dataset)
ds <- DESeqDataSetFromMatrix(countData=data, colData=meta_data, design=~meta_classes)

#Perform DEseq and create results object
Ds <- DESeq(ds)
res <- results(Ds)
print(res)

#Output results in.txt format
write.table(res,'DESeq2_PDX_R.txt', row.names=TRUE,col.names=NA, sep='\t', quote=TRUE, append=TRUE)

#Create MA plot
plotMA(res)

#Create enhance volcano plot with identified differentially expressed genes
EnhancedVolcano(res, x = 'log2FoldChange', y = 'pvalue', lab = rownames(res), xlim = c(-5, 8))
