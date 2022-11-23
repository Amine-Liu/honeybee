######################
####DESeq analysis####
######################

rm(list = ls())
options(stringsAsFactors = FALSE)
NE1 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/NE1.count", sep="\t", col.names = c("gene_name","NE1"))
NE2 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/NE2.count", sep="\t", col.names = c("gene_name","NE2"))
NE3 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/NE3.count", sep="\t", col.names = c("gene_name","NE3"))
N1 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/N1.count", sep="\t", col.names = c("gene_name","N1"))
N2 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/N2.count", sep="\t", col.names = c("gene_name","N2"))
N3 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/N2.count", sep="\t", col.names = c("gene_name","N3"))
F1 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/F1.count", sep="\t", col.names = c("gene_name","F1"))
F2 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/F2.count", sep="\t", col.names = c("gene_name","F2"))
F3 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/F3.count", sep="\t", col.names = c("gene_name","F3"))
D1 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/D1.count", sep="\t", col.names = c("gene_name","D1"))
D2 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/D2.count", sep="\t", col.names = c("gene_name","D2"))
D3 <- read.table("mypath/Apis_RNA_seq/output/htseq_count/count/D3.count", sep="\t", col.names = c("gene_name","D3"))

Worker=cbind(NE1,NE2[,2],NE3[,2],N1[,2],N2[,2],N3[,2],F1[,2],F2[,2],F3[,2],D1[,2],D2[,2],D3[,2])
tail(Worker)
Worker_filter =Worker[c(-1,-12321:-12325),]
tail(Worker_filter)
names(Worker_filter) <-c("gene_name","NE1","NE2","NE3","N1","N2","N3","F1","F2","F3","D1","D2","D3")
head(Worker_filter)
A<-Worker_filter
write.csv(A, "mypath/project/Apis_OSN_bulk_RNA-seq/_MT/rawcount_Worker_filter.csv")

library(DESeq2)
setwd("mypath/project/Apis_OSN_bulk_RNA-seq/_MT/")
rm(list = ls())
mycounts <- read.csv("mypath/project/Apis_OSN_bulk_RNA-seq/_MT/rawcount_Worker_filter.csv", header = T, row.names = 1)
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
condition <- factor(c(rep("NE", 3), rep("N", 3), rep("F", 3), rep("D", 3)), levels = c("NE","N","F","D"))
condition
colData <- data.frame(row.names = colnames(mycounts), condition)
colData
dds1 <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)
nrow(dds1)
dds1 <- dds1[rowSums(counts(dds1)) > 10, ]
nrow(dds1)
dds1 <- DESeq(dds1)


normalized_count<-counts(dds,normalized=TRUE)
write.csv(normalized_count,file = "mypath/Apis_RNA_seq/_MT/DESeq_normalized_DEG_worker.csv")
resultsNames(dds)
res = results(dds)
res = res[order(res$padj),]
head(res)
tail(res)
summary(res)
res<- as.data.frame(res)
res <- na.omit(res)
table(res$padj<0.01)

res = results(dds, contrast=c("condition", "N","NE"))
up_N_vs_NE_DEG <- subset(res, padj < 0.01 & log2FoldChange > 1.5)
dw_N_vs_NE_DEG <- subset(res, padj < 0.01 & log2FoldChange < -1.5)
write.csv(res, file="mypath/Apis_RNA_seq/_MT/N_vs_NE_DEG.csv")
write.csv(up_N_vs_NE_DEG, file="mypath/Apis_RNA_seq/_MT/up_N_vs_NE_DEG.csv")
write.csv(dw_N_vs_NE_DEG, file="mypath/Apis_RNA_seq/_MT/dw_N_vs_NE_DEG.csv")

res = results(dds, contrast=c("condition", "F","NE"))
up_F_vs_NE_DEG <- subset(res, padj < 0.01 & log2FoldChange > 1.5)
dw_F_vs_NE_DEG <- subset(res, padj < 0.01 & log2FoldChange < -1.5)
write.csv(res, file="mypath/Apis_RNA_seq/_MT/F_vs_NE_DEG.csv")
write.csv(up_F_vs_NE_DEG, file="mypath/Apis_RNA_seq/_MT/up_F_vs_NE_DEG.csv")
write.csv(dw_F_vs_NE_DEG, file="mypath/Apis_RNA_seq/_MT/dw_F_vs_NE_DEG.csv")

res = results(dds, contrast=c("condition", "D","NE"))
up_D_vs_NE_DEG <- subset(res, padj < 0.01 & log2FoldChange > 1.5)
dw_D_vs_NE_DEG <- subset(res, padj < 0.01 & log2FoldChange < -1.5)
write.csv(res, file="mypath/Apis_RNA_seq/_MT/D_vs_NE_DEG.csv")
write.csv(up_D_vs_NE_DEG, file="mypath/Apis_RNA_seq/_MT/up_D_vs_NE_DEG.csv")
write.csv(dw_D_vs_NE_DEG, file="mypath/Apis_RNA_seq/_MT/dw_D_vs_NE_DEG.csv")

res = results(dds, contrast=c("condition", "D","N"))
up_D_vs_N_DEG <- subset(res, padj < 0.01 & log2FoldChange > 1.5)
dw_D_vs_N_DEG <- subset(res, padj < 0.01 & log2FoldChange < -1.5)
write.csv(res, file="mypath/Apis_RNA_seq/_MT/D_vs_N_DEG.csv")
write.csv(up_D_vs_N_DEG, file="mypath/Apis_RNA_seq/_MT/up_D_vs_N_DEG.csv")
write.csv(dw_D_vs_N_DEG, file="mypath/Apis_RNA_seq/_MT/dw_D_vs_N_DEG.csv")

res = results(dds, contrast=c("condition", "D","F"))
up_D_vs_F_DEG <- subset(res, padj < 0.01 & log2FoldChange > 1.5)
dw_D_vs_F_DEG <- subset(res, padj < 0.01 & log2FoldChange < -1.5)
write.csv(res, file="mypath/Apis_RNA_seq/_MT/D_vs_F_DEG.csv")
write.csv(up_D_vs_F_DEG, file="mypath/Apis_RNA_seq/_MT/up_D_vs_F_DEG.csv")
write.csv(dw_D_vs_F_DEG, file="mypath/Apis_RNA_seq/_MT/dw_D_vs_F_DEG.csv")

res = results(dds, contrast=c("condition", "F","N"))
up_F_vs_N_DEG <- subset(res, padj < 0.01 & log2FoldChange > 1.5)
dw_F_vs_N_DEG <- subset(res, padj < 0.01 & log2FoldChange < -1.5)
write.csv(res, file="mypath/Apis_RNA_seq/_MT/F_vs_N_DEG.csv")
write.csv(up_F_vs_N_DEG, file="mypath/Apis_RNA_seq/_MT/up_F_vs_N_DEG.csv")
write.csv(dw_F_vs_N_DEG, file="mypath/Apis_RNA_seq/_MT/dw_F_vs_N_DEG.csv")

#merge DEG
DEG1<-read.csv("mypath/Apis_RNA_seq/_MT/N_vs_NE_DEG.csv", header = T, row.names = 1)
DEG2<-read.csv("mypath/Apis_RNA_seq/_MT/F_vs_NE_DEG.csv", header = T, row.names = 1)
DEG3<-read.csv("mypath/Apis_RNA_seq/_MT/D_vs_NE_DEG.csv", header = T, row.names = 1)
DEG4<-read.csv("mypath/Apis_RNA_seq/_MT/D_vs_N_DEG.csv", header = T, row.names = 1)
DEG5<-read.csv("mypath/Apis_RNA_seq/_MT/D_vs_F_DEG.csv", header = T, row.names = 1)
DEG6<-read.csv("mypath/Apis_RNA_seq/_MT/F_vs_N_DEG.csv", header = T, row.names = 1)
A<-rbind(DEG1, DEG2, DEG3, DEG4, DEG5, DEG6)
head(A)
diff_gene<-rownames(A)
diff_gene<-union(rownames(A),rownames(A))
head(diff_gene)
B<-read.csv("mypath/Apis_RNA_seq/_MT/DESeq_normalized_Q_W.csv",header = T, row.names = 1)
head(B)
head(match(diff_gene,rownames(B)))
C<-B[match(diff_gene,rownames(B)),]
head(C)
C<-data.matrix(C)
C<-na.omit(C)
dim(C)
head(C)
C<-C[,-(1:2)]
write.csv(C, file="mypath/Apis_RNA_seq/_MT/DESeq_normalized_DEG_worker.csv")

##################
####PCA_Worker####
##################

library(ggplot2)
Worker <- vst(dds1,blind=FALSE)
Worker
pdf("mypath/project/Apis_OSN_bulk_RNA-seq/_MT/Figures/DESeq_Worker_heatmap_PCA.pdf",width=5, height=4)
pcaData <- plotPCA(Worker, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()+
  theme(axis.line = element_line(colour="black"),plot.title = element_text(hjust = 0.5,size=16),
    panel.grid.minor = element_blank(),panel.background = element_rect(0, linetype = 0))
dev.off()

################################
####sample-to-sample heatmap####
################################

library(pheatmap)
library("RColorBrewer")
pdf("mypath/project/Apis_OSN_bulk_RNA-seq/_MT/Figures/sample-to-sample.pdf",width=6, height=5)
sampleDists <- dist(t(assay(Worker)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(Worker$condition, Worker$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,cluster_cols=FALSE,cluster_rows=FALSE,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

###############
####Volcano####
###############

rm(list = ls())
library(ggplot2)
pdf(file="~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/worker_volcano1.pdf")
F_vs_N <- read.csv(file = "~/project/Apis_OSN_bulk_RNA-seq/_MT/F_vs_N_DEG.csv", header = TRUE, sep =",")
D_vs_N <- read.csv(file = "~/project/Apis_OSN_bulk_RNA-seq/_MT/D_vs_N_DEG.csv", header = TRUE, sep =",")
D_vs_F <- read.csv(file = "~/project/Apis_OSN_bulk_RNA-seq/_MT/D_vs_F_DEG.csv", header = TRUE, sep =",")
NFD_vs_NE <- read.csv(file = "~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/NFD_vs_NE.csv", header = TRUE, sep =",")
cut_off_padj = 0.01
cut_off_log2FoldChange = 1.5

F_vs_N$change = ifelse(F_vs_N$padj < cut_off_padj & abs(F_vs_N$log2FoldChange) >= cut_off_log2FoldChange,
                          ifelse(F_vs_N$log2FoldChange > cut_off_log2FoldChange ,'up-regulated','down-regulated'),
                          'not-siginficant')

D_vs_N$change = ifelse(D_vs_N$padj < cut_off_padj & abs(D_vs_N$log2FoldChange) >= cut_off_log2FoldChange,
                          ifelse(D_vs_N$log2FoldChange > cut_off_log2FoldChange ,'up-regulated','down-regulated'),
                          'not-siginficant')

D_vs_F$change = ifelse(D_vs_F$padj < cut_off_padj & abs(D_vs_F$log2FoldChange) >= cut_off_log2FoldChange,
                          ifelse(D_vs_F$log2FoldChange > cut_off_log2FoldChange ,'up-regulated','down-regulated'),
                          'not-siginficant')
NFD_vs_NE$change = ifelse(NFD_vs_NE$padj < cut_off_padj & abs(NFD_vs_NE$log2FoldChange) >= cut_off_log2FoldChange,
                          ifelse(NFD_vs_NE$log2FoldChange > cut_off_log2FoldChange ,'up-regulated','down-regulated'),
                          'not-siginficant')


ggplot(
  F_vs_N, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue","grey", "red"))+
  scale_y_continuous(limits = c(0,50))+
  scale_x_continuous(limits = c(-5,5))+
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_padj),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)") + ggtitle("F_vs_N")+theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank())

ggplot(
  D_vs_N, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue","grey", "red"))+
  scale_y_continuous(limits = c(0,50))+
  scale_x_continuous(limits = c(-5,5))+
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_padj),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)") + ggtitle("D_vs_N")+theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank())

ggplot(
  D_vs_F, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue","grey", "red"))+
  scale_y_continuous(limits = c(0,10))+
  scale_x_continuous(limits = c(-5,5))+
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_padj),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)") + ggtitle("D_vs_F")+theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank())

ggplot(
  NFD_vs_NE, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue","grey", "red"))+
  scale_y_continuous(limits = c(0,60))+
  scale_x_continuous(limits = c(-5,5))+
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_padj),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)") + ggtitle("NFD_vs_NE")+theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank())
dev.off()

########################
####heatmap-all-DEGs####
########################

library(pheatmap)
library("RColorBrewer")
library(ComplexHeatmap)
library(ggplot2)
library(Hmisc)
pdf(file="~/project/Apis_OSN_bulk_RNA-seq/_MT/Figures/Figure2.pdf")
all_DEGs=read.csv("~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/TPM_DEG1.5_worker.csv",head=T, row.names = 1)
all_DEGs=as.matrix(all_DEGs)
range(all_DEGs)
pheatmap(all_DEGs,main="all_DEGs",
           cluster_cols = F,scale="row",
           cellwidth = 30,
           border_color=NA,
           show_rownames = F,
           angle_col ="0",
           treeheight_row = 0,)

#############################
####R-P relationship plot####
#############################

DEG1.5_worker=read.csv("~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/TPM_DEG1.5_worker_spearman.csv")
for(i in c(2:1055)){
matrix=data.frame(group=t(DEG1.5_worker[1,2:13]),Inclevel=t(DEG1.5_worker[i,2:13]))
matrix=as.matrix(matrix)
test=rcorr(matrix,type="spearman")
DEG1.5_worker[i,14]=test[[1]][2]
DEG1.5_worker[i,15]=test[[3]][2]}

DEG1.5_worker$change = ifelse(DEG1.5_worker$P<0.01,'Development',
                               ifelse(abs(DEG1.5_worker$R) < 0.3,'Division','Undefined'))
ggplot(data=DEG1.5_worker,mapping=aes(x=R,y=-log10(P),color=change))+
  geom_point(alpha=0.4, size=3.5)+
  scale_color_manual(values=c("#A0D600FF","#FF6D00FF","#d2dae2"))+
  geom_vline(xintercept=c(-0.3,0.3),col="red",lwd=0.8,linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01),col="blue",lwd=0.8,linetype = "dashed") +
  labs(x="R",y="-log10 (P)")+
  theme_bw()+
  ggtitle("R-P")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", legend.title = element_blank(),panel.grid=element_blank())

########################
####heatmap-division####
########################

DEG1.5_worker<-DEG1.5_worker[-1,]
DEG1.5_worker=DEG1.5_worker[order(DEG1.5_worker$R),]
DEG1.5_worker$order=c(1:nrow(DEG1.5_worker))
division <- DEG1.5_worker[which(abs(DEG1.5_worker$R)< 0.3),2:13]
rownames(division) <- DEG1.5_worker[which(abs(DEG1.5_worker$R)< 0.3),1]
data <- data.frame(NE=rowMeans(division[,1:3]),N=rowMeans(division[,4:6]),F=rowMeans(division[,7:9]),D=rowMeans(division[,10:12]))
rownames(data) <- rownames(division)
data=as.matrix(data)
map_div=ComplexHeatmap::pheatmap(data, main="division",
           cluster_cols = F,scale="row",
           cellwidth = 30,show_rownames = F,
           angle_col ="0",
           treeheight_row = 0, )



markgene <- c("LOC413754","LOC409576","CPR11","Hex70a","Mrjp1","Mrjp3","Mrjp4","Mrjp5","Mrjp7","LOC411290","LOC113218767","LOC100577778")
markgene <- as.data.frame(markgene)
map_div + rowAnnotation(link = anno_mark(at = which(rownames(data) %in% markgene$markgene),
    labels =rownames(data)[which(rownames(data) %in% markgene$markgene)], labels_gp = gpar(fontsize = 10)))


###########################
####heatmap-development####
###########################

dev <- DEG1.5_worker[which(abs(DEG1.5_worker$R)>0.9),2:13]
rownames(dev) <- DEG1.5_worker[which(abs(DEG1.5_worker$R)>0.9),1]
data <- data.frame(NE=rowMeans(dev[,1:3]),N=rowMeans(dev[,4:6]),F=rowMeans(dev[,7:9]),D=rowMeans(dev[,10:12]))
rownames(data) <- rownames(dev)
dev=as.matrix(data)
map_dev=ComplexHeatmap::pheatmap(dev,main="Development",
           cluster_cols = F,scale="row",
           cellwidth = 30,show_rownames = F,
           angle_col ="0",
           treeheight_row = 0, )



markgene <- c("LOC410021","LOC410368","Mrjp9","LOC100576135","Grp","LOC727165")
markgene <- as.data.frame(markgene)
map_dev + rowAnnotation(link = anno_mark(at = which(rownames(dev) %in% markgene$markgene),
    labels =rownames(dev)[which(rownames(dev) %in% markgene$markgene)], labels_gp = gpar(fontsize = 10)))

dev.off()

######################
####number of DEGs####
######################

library(dplyr)
library(ggplot2)
df <- data.frame(DEGs= c("up","up","up","up","down","down","down","down"),division_of_labor=c('NFD_vs_NE',"F_vs_N",'D_vs_N','D_vs_F','NFD_vs_NE',"F_vs_N",'D_vs_N','D_vs_F'),counts = c(415,39,25,5,295,49,49,7))
#df$counts <- as.numeric(df$counts)
counts_df <- df %>%
  group_by(DEGs,division_of_labor) %>%
  summarise(counts=n()) %>%
  mutate(normalized_counts = counts / sum(counts))
  counts_df <- counts_df %>%
  group_by(division_of_labor) %>%
  mutate(Virus_percentage_each_cell_type = normalized_counts / sum(normalized_counts)*100)
df$division_of_labor <- factor(df$division_of_labor,levels=c('NFD_vs_NE',"F_vs_N",'D_vs_N','D_vs_F'))
library(RColorBrewer)
newpalette <- c(brewer.pal(9, "YlOrRd")[c(4)],brewer.pal(9, "Blues")[c(4)])
pdf("~/project/Apis_OSN_bulk_RNA-seq/_MT/Figures/number-cluster-DEG.pdf",width = 10)

ggplot(data=df,aes(x=division_of_labor,y=counts,fill=DEGs))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+
  scale_fill_manual(values=newpalette)+
  labs(x="division_of_labor",y="Number of DEGs",fill="DEGs",title="Number of DEGs") +
  theme(axis.line = element_line(colour="black"),axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5,size=16),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_rect(0, linetype = 0))
dev.off()


DEG1<-read.csv("~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/D_vs_N_DEG1.5.csv", header = T, row.names = 1)
DEG2<-read.csv("~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/D_vs_F_DEG1.5.csv", header = T, row.names = 1)
DEG3<-read.csv("~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/F_vs_N_DEG2.csv", header = T, row.names = 1)
A<-rbind(DEG1, DEG2, DEG3)
head(A)
diff_gene<-rownames(A)
diff_gene<-union(rownames(A),rownames(A))
head(diff_gene)
B<-read.csv("~/project/Apis_OSN_bulk_RNA-seq/_MT/TPM_Worker.csv",header = T, row.names = 1)
head(B)
head(match(diff_gene,rownames(B)))
C<-B[match(diff_gene,rownames(B)),]
head(C)
C<-data.matrix(C)
C<-na.omit(C)
dim(C)
head(C)
write.csv(C, file="~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/padj0.01_1.5_diffgeneNFD_list.csv")


####################
####qPCR barplot####
####################
library(ggplot2)
test <- read.csv("~/project/Apis_OSN_bulk_RNA-seq/_MT/bulkRNA-seq/DEGs.csv", header = T, row.names = 1)
test$division=c(rep("NE",3),rep("N",3),rep("F",3),rep("D",3))
a=list()

for (i in 1:12) {
mean<-aggregate(test[,i],by=list(test$division),FUN=mean)
sd<-aggregate(test[,i],by=list(test$division),FUN=sd)
N<-aggregate(test[,i],by=list(test$division),FUN=length)
name=colnames(test)[i]
test1=data.frame(mean,sd$x,N$x)
colnames(test1)=c("Division","mean","sd","N")
test1$se=test1$sd / sqrt(test1$N)
test1$Division=factor(test1$Division,levels=c("NE","N","F","D"))
a[[i]]=ggplot(test1, aes(x=Division, y=mean, fill=Division))+
      geom_bar(color="black",stat="identity",width=.6)+
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,position=position_dodge(.6))+
      theme_bw()+theme(panel.grid=element_blank())+
      ggtitle(name)+labs(y="Percent of inclusion form")
}


pdf("~/project/Apis_OSN_bulk_RNA-seq/_MT/Figures/qPCR_DEGs_hist.pdf")
a[[1]]
a[[2]]
a[[3]]
a[[4]]
a[[5]]
a[[6]]
a[[7]]
a[[8]]
a[[9]]
a[[10]]
a[[11]]
a[[12]]

dev.off()
