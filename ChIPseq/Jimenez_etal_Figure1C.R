#Script used to generate plot seen in Figure 1C
#Files available on GEO GSE169250 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169250
#Convert bigwig files to bedgraph prior to running files through this script

#Load ggplot2 for plotting
library(ggplot2)

#Read in ChIP bedgraph files and corresponding peak bed files
LS01_chip<-read.table("DPY-27_LS01_emb_avg_LAS06_LAS07_ratio.bdg")
COP325_chip<-read.table("DPY-27_COP325_emb_ext120_ext122_avg_ratio.bdg")
LS01_peaks<-read.table("DPY-27_LS01_emb_avg_LAS06_LAS07_chip_peaks.bed")
COP325_peaks<-read.table("DPY-27_COP325_emb_avg_LW15_LW24_LW35_chip_peaks.bed")

#Subset peaks on chromosomes of interest to new tables
LS01chrIIpeaks<-LS01_peaks[which(LS01_peaks$V1=='chrII'),]
LS01chrXpeaks<-LS01_peaks[which(LS01_peaks$V1=='chrX'),]
COP325chrIpeaks<-COP325_peaks[which(COP325_peaks$V1=='chrI'),]
COP325chrXpeaks<-COP325_peaks[which(COP325_peaks$V1=='chrX'),]

#Define regions of interest for each strain & chromosome. In this case a 200kb region surrounding the centerpoint of insertion
COP325_suba<-4660000-100000
COP325_subb<-4660000+100000
LS01_suba<-8420000-100000
LS01_subb<-8420000+100000

#Assign strain names to each ChIP dataset prior to creating a master dataset
LS01_chip$strain<-"LS01"
COP325_chip$strain<-"COP325"

#Defines a function to assign peaks to a dataset. If a datapoint lies within the coordinates of a defined peak, it is classified as a peak in the dataframe
assign.peaks<-function(dataset,peaksfile){
  dataset$group<-"No"
  for(i in 1:length(peaksfile$V2)){
    dataset$group[which(dataset$V1==peaksfile$V1[i] & dataset$V3>peaksfile$V2[i] & dataset$V3<peaksfile$V3[i])]<-"Yes"
  }
  output<-dataset
  return(output)
}

#Runs function on each dataset
LS01_chip_peaks<-assign.peaks(LS01_chip,LS01_peaks)
COP325_chip_peaks<-assign.peaks(COP325_chip,COP325_peaks)

#Combine ChIP datasets into one master dataframe
Master_Chip<-rbind(LS01_chip_peaks,COP325_chip_peaks)

#Removes mitochondrial DNA
Master_Chip<-Master_Chip[which(Master_Chip$V1!="chrM"),]

#Subsets ChrX from master dataframe and defines class for plotting
Sub_chipa<-Master_Chip[which(Master_Chip$V1=="chrX"),]
Sub_chipa$class<-"ChrX"

#Subsets ChrII data from LS01 strain data that lies within user defined region of interest from master dataframe and defines class for plotting
Sub_chipb<-Master_Chip[which(Master_Chip$strain=="LS01" & Master_Chip$V1=='chrII' & Master_Chip$V2>LS01_suba & Master_Chip$V3<LS01_subb),]
Sub_chipb$class<-"ChrII 200Kb"

#Subsets ChrI data from COP325 strain data that lies within user defined region of interest from master dataframe and defines class for plotting
Sub_chipd<-Master_Chip[which(Master_Chip$strain=="COP325" & Master_Chip$V1=='chrI' & Master_Chip$V2>COP325_suba & Master_Chip$V3<COP325_subb),]
Sub_chipd$class<-"ChrI 200Kb"

#Combine filtered data into one dataframe
Sub_chip<-rbind(Sub_chipa,Sub_chipb,Sub_chipd)

#Script to generate plot in Figure 1C
ggplot(Sub_chip,aes(x=strain,y=V4,color=class))+
  geom_violin()+
  coord_cartesian(ylim = c(0,3))
