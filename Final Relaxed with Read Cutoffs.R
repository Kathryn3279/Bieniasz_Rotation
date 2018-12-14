
#Read in RNAseq data
RNAseq_data1<-read.csv("C:/Users/kathr/Desktop/Rotations/Bieniasz Rotation/rowCount_Matrix.csv")
#Subset averages for each gene
RNAseq_Ave<-RNAseq_data1[,c(1,26:31)]

#Calculate each condition's reads as a percent of total reads
RNAseq_Ave$RSVpercents<-(RNAseq_Ave$RSV*100)/(RNAseq_Ave$RSV+RNAseq_Ave$RSV.GAMMA+RNAseq_Ave$MOCK.GAMMA+RNAseq_Ave$RSV.ALPHA+RNAseq_Ave$MOCK.ALPHA)
RNAseq_Ave$RSVGammapercents<-(RNAseq_Ave$RSV.GAMMA*100)/(RNAseq_Ave$RSV+RNAseq_Ave$RSV.GAMMA+RNAseq_Ave$MOCK.GAMMA+RNAseq_Ave$RSV.ALPHA+RNAseq_Ave$MOCK.ALPHA)
RNAseq_Ave$MGammapercents<-(RNAseq_Ave$MOCK.GAMMA*100)/(RNAseq_Ave$RSV+RNAseq_Ave$RSV.GAMMA+RNAseq_Ave$MOCK.GAMMA+RNAseq_Ave$RSV.ALPHA+RNAseq_Ave$MOCK.ALPHA)
RNAseq_Ave$MAlphapercents<-(RNAseq_Ave$MOCK.ALPHA*100)/(RNAseq_Ave$RSV+RNAseq_Ave$RSV.GAMMA+RNAseq_Ave$MOCK.GAMMA+RNAseq_Ave$RSV.ALPHA+RNAseq_Ave$MOCK.ALPHA)
RNAseq_Ave$RSVAlphapercents<-(RNAseq_Ave$RSV.ALPHA*100)/(RNAseq_Ave$RSV+RNAseq_Ave$RSV.GAMMA+RNAseq_Ave$MOCK.GAMMA+RNAseq_Ave$RSV.ALPHA+RNAseq_Ave$MOCK.ALPHA)

CombinedList2<-list()

#Looping through several percent cutoffs(x)
x<-seq(50,85,by=5)
i=50
for (i in x){
  #Finding genes with average reads in any condition great than or equal to x percent//"uniques"
  AllURSV<-as.vector(which(RNAseq_Ave$RSVpercents >= i ,arr.ind = TRUE))
  AllURSVGamma<-as.vector(which(RNAseq_Ave$RSVGammapercents >= i ,arr.ind = TRUE))
  AllUMGamma<-as.vector(which(RNAseq_Ave$MGammapercents >= i ,arr.ind = TRUE))
  AllUMAlpha<-as.vector(which(RNAseq_Ave$MAlphapercents >= i ,arr.ind = TRUE))
  AllURSVAlpha<-as.vector(which(RNAseq_Ave$RSVAlphapercents >= i ,arr.ind = TRUE))
  
  library(tidyverse)
  y<-((100-i)/5)-1
  z<-(100-i)
  BadURSV<-as.vector(which(between(RNAseq_Ave$RSVpercents,y,z),arr.ind = TRUE))
  BadURSVGamma<-as.vector(which(between(RNAseq_Ave$RSVGammapercents,y,z),arr.ind = TRUE))
  BadUMGamma<-as.vector(which(between(RNAseq_Ave$MGammapercents,y,z),arr.ind = TRUE))
  BadUMAlpha<-as.vector(which(between(RNAseq_Ave$MAlphapercents,y,z),arr.ind = TRUE))
  BadURSVAlpha<-as.vector(which(between(RNAseq_Ave$RSVAlphapercents,y,z),arr.ind = TRUE))
  
  BadAll<-as.vector(unique(c(BadURSV,BadURSVGamma,BadUMGamma,BadUMAlpha,BadURSVAlpha)))
  
  #Find genes with condition ave reads greater than 50 in sum not average
  RNAseq_data1$RSVsum<-RNAseq_data1$Sample_1B7_RSV+RNAseq_data1$Sample_1C10_RSV+RNAseq_data1$Sample_2E10_RSV+RNAseq_data1$Sample_3G7_RSV
  RNAseq_data1$RSVGammasum<-RNAseq_data1$Sample_1B7_RSV_GAMMA+RNAseq_data1$Sample_1C10_RSV_GAMMA+RNAseq_data1$Sample_2E10_RSV_GAMMA+RNAseq_data1$Sample_3G7_RSV_GAMMA
  RNAseq_data1$MockGammasum<-RNAseq_data1$Sample_1B7_MOCK_GAMMA+RNAseq_data1$Sample_1C10_MOCK_GAMMA+RNAseq_data1$Sample_2E10_MOCK_GAMMA+RNAseq_data1$Sample_3G7_MOCK_GAMMA
  RNAseq_data1$MockAlphasum<-RNAseq_data1$Sample_1B7_MOCK_ALPHA+RNAseq_data1$Sample_1C10_MOCK_ALPHA+RNAseq_data1$Sample_2E10_MOCK_ALPHA+RNAseq_data1$Sample_3G7_MOCK_ALPHA
  RNAseq_data1$RSVAlphasum<-RNAseq_data1$Sample_1B7_RSV_ALPHA+RNAseq_data1$Sample_1C10_RSV_ALPHA+RNAseq_data1$Sample_2E10_RSV_ALPHA+RNAseq_data1$Sample_3G7_RSV_ALPHA
  
  
  RSVhighReads<-as.vector(which(RNAseq_data1$RSVsum > 50, arr.ind=TRUE))
  RSVGammahighReads<-as.vector(which(RNAseq_data1$RSVGammasum > 50, arr.ind=TRUE))
  MGammahighReads<-as.vector(which(RNAseq_data1$MockGammasum > 50, arr.ind=TRUE))
  MAlphahighReads<-as.vector(which(RNAseq_data1$MockAlphasum > 50, arr.ind=TRUE))
  RSVAlphahighReads<-as.vector(which(RNAseq_data1$RSVAlphasum > 50, arr.ind=TRUE))
  
  #find intersect
  URSVtoKeep<-intersect(RSVhighReads,AllURSV)
  URSVGammatoKeep<-intersect(RSVGammahighReads,AllURSVGamma)
  UMGammatoKeep<-intersect(MGammahighReads,AllUMGamma)
  UMAlphatoKeep<-intersect(MAlphahighReads,AllUMAlpha)
  URSVAlphatoKeep<-intersect(RSVAlphahighReads,AllURSVAlpha)
  
  #Find gene name for uniques
  UniqueRSVGene<-as.vector(RNAseq_Ave[URSVtoKeep,1]) 
  UniqueRSVGammaGene<-as.vector(RNAseq_Ave[URSVGammatoKeep,1]) 
  UniqueMGammaGene<-as.vector(RNAseq_Ave[UMGammatoKeep,1])
  UniqueMAlphaGene<-as.vector(RNAseq_Ave[UMAlphatoKeep,1]) 
  UniqueRSVAlphaGene<-as.vector(RNAseq_Ave[URSVAlphatoKeep,1]) 
  
  #find gene name for "Bad uniques"
  BadAllGenes<-as.vector(RNAseq_Ave[BadAll,1]) 
  
  #Write table for final unique gene names and cutoffs
  FinalRSVList<-as.list(setdiff(UniqueRSVGene,BadAllGenes))
  FinalRSVList<-append(FinalRSVList,i,after=0)
  FinalRSVList<-append(FinalRSVList,"RSV",after=0)

  
  FinalRSVGammaList<-as.list(setdiff(UniqueRSVGammaGene,BadAllGenes))
  FinalRSVGammaList<-append(FinalRSVGammaList,i,after=0)
  FinalRSVGammaList<-append(FinalRSVGammaList,"RSVGamma",after=1)

  
  FinalMGammaList<-as.list(setdiff(UniqueMGammaGene,BadAllGenes))
  FinalMGammaList<-append(FinalMGammaList,"MockGamma",after=0)
  FinalMGammaList<-append(FinalMGammaList,i,after=1)
  
  FinalMAlphaList<-as.list(setdiff(UniqueMAlphaGene,BadAllGenes))
  FinalMAlphaList<-append(FinalMAlphaList,i,after=0)
  FinalMAlphaList<-append(FinalMAlphaList,"MockAlpha",after=1)

  
  FinalRSVAlphaList<-as.list(setdiff(UniqueRSVAlphaGene,BadAllGenes))
  FinalRSVAlphaList<-append(FinalRSVAlphaList,i,after=0)
  FinalRSVAlphaList<-append(FinalRSVAlphaList,"RSVAlpha",after=1)
  
  FullList<-as.character(c(FinalMAlphaList,FinalMGammaList,FinalRSVList,FinalRSVAlphaList,FinalRSVGammaList))
  
  CombinedList2<-append(CombinedList2,FullList,after=length(CombinedList2))
    
}
#Writing results to tables
CombinedList2<-as.character(CombinedList2)
write.csv(CombinedList2,file="TestCombined.csv")
