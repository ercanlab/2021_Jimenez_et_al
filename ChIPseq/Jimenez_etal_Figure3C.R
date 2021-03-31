#Script used to generate plot seen in Figure 1C
#Files available on GEO GSE169250 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169250
#Convert bigwig files to bedgraph prior to running files through this script

#Load ggplot2 for plotting
library(ggplot2)

#Read in ChIP bedgraph files
#Read in rex site annotation file
N2_chip<-read.table("DPY27_N2_emb_avg_CJ132_CJ39_SE172_ratio.bdg")
LS05_chip<-read.table("DPY-27_LS05_emb_avg_LAS55_DJ84_ratio.bdg")
BR01_chip<-read.table("DPY27_BR01A_emb_avg_JK12_JK23_JK24_ratio.bdg")
rex_sites<-read.table("~/Desktop/Data/bed/rex1bp.bed")

#Remove unnecessary columns from rex annotation dataframe
rex_sites<-rex_sites[,c(1:5)]

#Define insertion site on chrX
chrX_ins<-as.data.frame(cbind("chrX",14373128,14373129,"chrX_insert","chrX_insert"))

#Make a loci dataframe containing rex sites and insertion site
loci<-as.data.frame(rbind(rex_sites,chrX_ins))

#Classify rex sites by distance from the insertion site
loci$group<-"NA"
loci$group[which(loci$V2>12373128 & loci$V2<16373128)]<-"4Mb"
loci$group[which(loci$V2>13373128 & loci$V2<15373128)]<-"2Mb"
loci$group[65]<-"Insert"

#Change data class for analysis 
loci$V2<-as.numeric(loci$V2)
loci$V3<-as.numeric(loci$V3)

#Define regions around each rex to consider for analysis
start1<-as.numeric(loci$V2)-250
stop1<-as.numeric(loci$V2)-100
start2<-as.numeric(loci$V2)+100
stop2<-as.numeric(loci$V2)+250

#Add start and stop for regions of interest to the loci dataframe
loci$start1<-start1
loci$start2<-start2
loci$stop1<-stop1
loci$stop2<-stop2

#Assign each dataframe strain information and combine dataframes into a master dataframe. Classify appropriate columns accordingly
N2_chip$strain<-"N2"
BR01_chip$strain<-"BR01"
LS05_chip$strain<-"LS05"
Master_Dataset<-as.data.frame(rbind(N2_chip,BR01_chip,LS05_chip))
Master_Dataset$V2<-as.numeric(Master_Dataset$V2)
Master_Dataset$V3<-as.numeric(Master_Dataset$V3)
Master_Dataset$V4<-as.numeric(Master_Dataset$V4)

#Define function to assign corresponding rex information to the data
assign.rex<-function(dataset,selection){
  dataset$rex<-as.character("NA")
  dataset$group<-as.character("NA")
  for(i in 1:length(selection$V2)){
    dataset$rex[which(dataset$V1==selection$V1[i] & dataset$V2>selection$start1[i] & dataset$V2<selection$stop2[i])]<-as.character(selection$V4[i])
    dataset$group[which(dataset$V1==selection$V1[i] & dataset$V2>selection$start1[i] & dataset$V2<selection$stop2[i])]<-as.character(selection$group[i])
  }
  output<-dataset
  return(output)
}

#Subset data which is on the X chromosome and run the above defined scripts to determine which data points to keep and assign rex information to master dataframe
Master_X_Dataset<-Master_Dataset[which(Master_Dataset$V1=="chrX"),]
Master_X_Dataset$V1<- factor(Master_X_Dataset$V1)
Master_X_Dataset_rex<-assign.rex(Master_X_Dataset,loci)

#Removes datapoints that were not assigned rex information
Master_X_Dataset_rex_final<-Master_X_Dataset_rex[which(Master_X_Dataset_rex$rex!="NA"),]

#Define function to assign which data points to filter out from analysis
filter_out<-function(dataset,selection){
  dataset$keep<-"Yes"
  for(i in 1:length(selection$V2)){
    dataset$keep[which(dataset$V1==selection$V1[i] & dataset$V2>selection$stop1[i] & dataset$V2<selection$start2[i])]<-"No"
  }
  output<-dataset
  return(output)
}

#Run filter_out function on data and remove data from analysis
Master_X_Dataset_rex_final_pre_filter<-filter_out(Master_X_Dataset_rex_final,loci)
Master_X_Dataset_rex_final_v2<-Master_X_Dataset_rex_final_pre_filter[which(Master_X_Dataset_rex_final_pre_filter$keep=="Yes"),]

#Subset data for processing
Final_WT<-Master_X_Dataset_rex_final_v2[which(Master_X_Dataset_rex_final_v2$strain=="N2"),]
Final_LS05<-Master_X_Dataset_rex_final_v2[which(Master_X_Dataset_rex_final_v2$strain=="LS05"),]
Final_BR01<-Master_X_Dataset_rex_final_v2[which(Master_X_Dataset_rex_final_v2$strain=="BR01"),]

#Calculate mean DCC binding level for each assigned rex for each strain
sub_WT<-0
sub_BR01<-0
sub_LS05<-0
for(i in 1:length(loci$V1)){
  sub_WT[i]<-mean(Final_WT$V4[which(Final_WT$rex==loci$V4[i])])
  sub_BR01[i]<-mean(Final_BR01$V4[which(Final_BR01$rex==loci$V4[i])])
  sub_LS05[i]<-mean(Final_LS05$V4[which(Final_LS05$rex==loci$V4[i])])
}

all_loci<-cbind(loci,sub_WT,sub_BR01,sub_LS05)

#Load ggplot2 for plotting
library(ggplot2)

#Plot correlation plots for both insertions
ggplot(all_loci,aes(x=sub_BR01,y=sub_LS05,color=group))+
  geom_point()+
  scale_color_manual(values=c("#1700FF","#1000B2","#FE2D00","#000000"))

#Plot correlation plots for downstream oriented insertion and wildtype
ggplot(all_loci,aes(x=sub_BR01,y=sub_WT,color=group))+
  geom_point()+
  scale_color_manual(values=c("#1700FF","#1000B2","#FE2D00","#000000"))

#Plot correlation plots for upstream oriented insertion and wildtype
ggplot(all_loci,aes(x=sub_LS05,y=sub_WT,color=group))+
  geom_point()+
  scale_color_manual(values=c("#1700FF","#1000B2","#FE2D00","#000000"))
