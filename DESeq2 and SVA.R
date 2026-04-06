### IBD ###
library(DESeq2)

Counts=read.table("E:/biodata/AH1.count")[,2]
for (i in dir("E:/biodata/",full.names = T)[-1]){
  Counts=cbind(Counts,read.table(i)[,2])
}
rownames(Counts)=read.table("E:/biodata/AH1.count")[,1]
colnames(Counts)=gsub(".count","",dir("E:/biodata/"))
View(Counts)
Counts=Counts[-grep("__",rownames(Counts)),]
write.csv(Counts, "E:/biodata/counts.csv")

library(sva)
mycount=read.table("E:/biodata/counts.txt", header = TRUE)
rownames(mycount)=mycount$Column1
mycount=mycount[,-1]
mycount=as.matrix(mycount)
batch <- c(rep(1, 14), rep(2, 53))
group <- c(rep(1, 14), rep(2, 53))
Crohns_Disease_normalized <- ComBat_seq(mycount, batch=batch, group=group, full_mod=FALSE)
write.table(Crohns_Disease_normalized, "E:/biodata/Crohns_Disease_normalized.txt")

ascending_colon_CD <- read.table("E:/biodata/ascending colon CD.txt", header = TRUE)
rownames(ascending_colon_CD)=ascending_colon_CD$Column1
ascending_colon_CD=ascending_colon_CD[,-1]
CountsMatrix=as.matrix(Counts)
View(CountsMatrix)
Class=read.csv("E:/biodata/GroupFiles.csv")[,2]
Class=data.frame(patients=as.factor(Class))
View(Class)
DESeq2Obj=DESeqDataSetFromMatrix(CountsMatrix,Class,~patients)
CountsNorm=DESeq(DESeq2Obj)
DEGresult=results(CountsNorm)
View(DEGresult)
CountsNorm$sizeFactor
DEGresult=data.frame(results(CountsNorm))
View(DEGresult)
DEGresult$padj=p.adjust(DEGresult$pvalue,method = "BH")
DEGresult=DEGresult[order(DEGresult$padj),]
View(DEGresult)
write.csv(DEGresult, "d:/New folder/MS/DEGS.csv")
Upregulated=DEGresult[which(DEGresult$log2FoldChange > 1 & DEGresult$pvalue < 0.05),]
View(UPDEG)
dim(Upregulated)
rownames(UPDEG)
write.csv(Upregulated,"E:/biodata/Upregulated.csv")
Downregulated=DEGresult[which(DEGresult$log2FoldChange < 1 & DEGresult$pvalue < 0.05),]
View(DOWNDEG)
dim(Downregulated)
rownames(DOWNDEG)
write.csv(Downregulated,"E:/biodata/Downregulated.csv")
