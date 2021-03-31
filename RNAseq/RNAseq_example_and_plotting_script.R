#Script to show how to use the functions defined in the RNAseq_data_processing_script

#Load ggplot2 for plotting, and stringr and dplyr for data processing steps
library(ggplot2)
library('stringr')
library(dplyr)

#To test the usage of the script and to run the plots, use the sample file provided in the github to generate the plots seen in the figures as outlined below
#Read in DESEQ output file & run the defined functions
LS01emb<-read.table("controlvstreatment.deseq.txt",header=T)
LS01emb<-get.chr(LS01emb)
LS01emb<-get.pos(LS01emb)

#At this point the data is ready for basic plotting. If you would like to combine datasets into one dataframe or further group your data to run statistical tests, now is where you can define the various groups.
LS01emb$strain<-"LS01emb"

#Grouping by chromosomes
LS01emb$XvA<-"A"
LS01emb$XvA[which(LS01emb$chr.name=='X')]<-"X"
LS01emb$XvA[which(LS01emb$chr.name=='II')]<-"II"
LS01emb$XvA[which(LS01emb$chr.name=='I')]<-"I"

#Using ggplot to generate the plot seen in Figure 1E.
pq <- ggplot(LS01emb, aes(x=strain, y=log2FoldChange, fill = XvA)) +
  geom_boxplot(linetype='dashed',outlier.shape = NA,notch=TRUE)+scale_fill_brewer(palette="Dark2") +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  coord_cartesian(ylim=c(-1.8,1.8))+
  ylab("log2FoldChange") + xlab("Chromosome") #+ geom_boxplot(outlier.shape =NA)
pq+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Statistical testing to determine statistical significance of log2fc comparing chrII, chrIV & chrV to chrI, chrII or chrX. 
t.test(LS01emb$log2FoldChange[which(LS01emb$XvA=="I")],LS01emb$log2FoldChange[which(LS01emb$XvA=="A")])
t.test(LS01emb$log2FoldChange[which(LS01emb$XvA=="II")],LS01emb$log2FoldChange[which(LS01emb$XvA=="A")])
t.test(LS01emb$log2FoldChange[which(LS01emb$XvA=="X")],LS01emb$log2FoldChange[which(LS01emb$XvA=="A")])

#Subsetting data and preparing data for gene by gene log2fc analysis.
LS01II<-LS01emb[which(LS01emb$chr.name=='II'),]
LS01II$chr.start<-as.numeric(as.character(LS01II$chr.start))
LS01II$chr.stop<-as.numeric(as.character(LS01II$chr.stop))

#Using ggplot to generate the plot seen in Figure 1D.
genebygene <- ggplot(LS01II, aes(x=(chr.start+chr.stop)/2, y=log2FoldChange)) +
  geom_point()+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE))+
  coord_cartesian(xlim=c(8420088-50000,8420088+49999))+
  geom_rect(aes(xmin=chr.start,ymin=-6,xmax=chr.stop,ymax=-5.5),alpha=.4)+
  xlab("Chromosome Coordinate")
genebygene 

#Subsetting data and preparing data for bin by bin log2fc analysis.
IIsub<-LS01emb[which(LS01emb$chr.name=="II"),]
IIsub$chr.start<-as.numeric(paste(IIsub$chr.start))

#Uses the central insertion as a centerpoint and assigns datapoints to bins sized according to width defined in the mutate line.
IIsub <- IIsub %>% 
  mutate( bin=cut_width(IIsub$chr.start, width=500000, center =84200000))

#Using ggplot to generate the plots seen in Figure S1B.
pq <- ggplot(IIsub, aes(x=bin, y=log2FoldChange)) +
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim=c(-2,2))+
  ylab("log2FoldChange") + xlab("Chromosome Position")
pq + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Defining which genes are upregulated and which genes are downregulated in each bin to run Fisher Exact test
all.up<-length(IIsub$log2FoldChange[which(IIsub$log2FoldChange>0)])
all.down<-length(IIsub$log2FoldChange[which(IIsub$log2FoldChange<0)])
sub.up<-c()
sub.down<-c()
pval<-c()
for(i in 1:length(levels(IIsub$bin))){
  sub.up[i]<-length(IIsub$log2FoldChange[which(IIsub$log2FoldChange>0 & IIsub$bin==levels(IIsub$bin)[i])])
  sub.down[i]<-length(IIsub$log2FoldChange[which(IIsub$log2FoldChange<0 & IIsub$bin==levels(IIsub$bin)[i])])
  test.fish<-matrix(c(sub.up[i],sub.down[i],all.up,all.down),ncol=2,byrow=T)
  pval[i]<-as.numeric(fisher.test(test.fish,alternative='two.sided')$p.value)
}
#Provides output of genes at user defined statistical significance level
levels(IIsub$bin)[which(pval<.001)]
