rm(list = ls(all.names = TRUE))
setwd("E:/Shirin_PourAkaber/counts") 
full_count <- read.table("counts_N.csv") 
# Remove Y chromosome genes with 0 expression 
full_count2 = full_count [-c(25,5510,5512,12088,12480,12515,12519, 12521,14519, 15182,15241, 15260, 15931, 15957,
                             16127, 17117, 17604, 17693, 20282 , 22914,24619 , 24726, 24749, 24792,24932,
                             26223, 26582, 27272, 28111, 28220, 28649, 29570, 32358, 32591, 33300, 33891, 
                             34011, 34350, 34536, 48501, 50946, 53498, 54421, 56690, 56844 ),]
head(full_count2)
# Check for duplicates
length(unique(rownames(full_count2))) == nrow(full_count2)
colnames(full_count2)=c("CVID49_N", "CVID50_N", "CVID46_N", "CVID47_N",
                        "CVID40_N", "CVID45_N", "CVID53_N", "CVID41_N",
                        "CVID51_N", "CVID52_N", "CVID48_N",
                        "HD1_N", "HD2_N", "HD3_N", "HD4_N")
full_count2[is.na(full_count2)] <- 0
# Bar plot of the count data
ctl_sum2=apply(full_count2,2,sum)
barplot(ctl_sum2,las=2)
boxplot(log(full_count2+1),las=2, par(cex.axis=0.6, cex.lab=1, cex.main=1.2, cex.sub=1))
# Differential expression analysis
library(DESeq2)
# Subtype= CVID with or without Autoimmune Deficiency
subtype=factor(c("CVID w/o. AD", "CVID w/o. AD", "CVID w. AD", "CVID w. AD",
                 "CVID w. AD", "CVID w. AD", "CVID w. AD", "CVID w. AD",
                 "CVID w/o. AD", "CVID w/o. AD", "CVID w/o. AD",
                 "HD","HD","HD","HD"),levels=c("CVID w/o. AD","CVID w. AD","HD"))
Disease2=factor(  c(rep("CVID",11),rep("HD",4))  ,levels=c("HD","CVID"))
CVID_HD_N=data.frame(sampleName=colnames(full_count2),Disease=Disease2, subtype=subtype)
# Design Formula
ddsMat.CVID_HD_N <- DESeqDataSetFromMatrix(countData = round(full_count2),colData = CVID_HD_N,design = ~ Disease)
dds.CVID_HD_N <- DESeq(ddsMat.CVID_HD_N)
# VST function
rld_CVID_HD_N<- vst(dds.CVID_HD_N) 
# See the results
res.CVID_HD_N <- results(dds.CVID_HD_N) # Extract results
head(res.CVID_HD_N)
summary(res.CVID_HD_N)
# Extract the normalized data and calculate the mean count of it (for identifing markers)
norm_cts.CVID_HD_N=counts(dds.CVID_HD_N, normalized=TRUE)
dim(res.CVID_HD_N)
# Adjusted p-value < 0.05 & remove NAs
DEGs_CVID_HD_N=res.CVID_HD_N[(res.CVID_HD_N$padj)<0.05 & !is.na(res.CVID_HD_N$padj),]
# PCA plot
m3= plotPCA(rld_CVID_HD_N[row.names(DEGs_CVID_HD_N)], intgroup=c("Disease","subtype"))
names(m3)
m3$data
# Add labels to the plots
library(ggplot2)
m3 <- m3 + geom_text(aes_string(x = m3$data[,1], y = m3$data[,2], label = "name"), color = "black",cex=2)
print(m3)
# Order the genes based on the log2FC
order_genes=DEGs_CVID_HD_N[order(abs(DEGs_CVID_HD_N$log2FoldChange),decreasing = TRUE),]
# Remove HLA genes
order_genes = order_genes [-c(3,4,7,8,12,14,16,17,22,43,57,66,70,74,75,80,81,85,88,93,108),]
head(order_genes)
dim(order_genes)
df=data.frame(order_genes)
# Collect the top 50 genes
top50=order_genes[1:50,]
head(top50)
names_top50=rownames(top50)

rlog_trans_cts=assay(rlog(dds.CVID_HD_N))
head(rlog_trans_cts)
counts_top50=rlog_trans_cts[names_top50,]
head(counts_top50)
tail(counts_top50)

library(pheatmap)
pheatmap(as.matrix(t(counts_top50)),scale="column",fontsize = 8,clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",clustering_method = "average")

# Pathway Analysis ( Transcription IDs were removed using a Macro in Excel)
library("biomaRt")
# Use Ensemble database
mart1<-  useMart("ensembl") 
mart<-useDataset("hsapiens_gene_ensembl",mart1)
values<-rownames(df)
biomart2=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","gene_biotype"), filters = "ensembl_gene_id", values = values, mart= mart)
# To convert the ENSEMBEL ids to ENTREZ Ids:
library(edgeR)
library(org.Hs.eg.db)
library(GO.db)
dfnames<-row.names(df)
class(dfnames)
annot=select(org.Hs.eg.db,keys=dfnames,columns=c("ENTREZID"), keytype="ENSEMBL")
entrezID=annot$ENTREZID
entrezID2 = entrezID[!is.na(entrezID)]
go <- goana(entrezID2,species="Hs")
dim(go)
topGO <-topGO(go,ontology="BP" ,n=729)
head(topGO)
keg <- kegga(entrezID2, species="Hs")
head(keg)
topKEGG <- topKEGG(keg)
