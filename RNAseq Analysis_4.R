RNAseq Analysis_4
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

#Finding genes with average reads in any condition great than or equal to 90 percent//"uniques"
AllURSV<-as.vector(which(RNAseq_Ave$RSVpercents >= 90.0 ,arr.ind = TRUE))
AllURSVGamma<-as.vector(which(RNAseq_Ave$RSVGammapercents >= 90.0 ,arr.ind = TRUE))
AllUMGamma<-as.vector(which(RNAseq_Ave$MGammapercents >= 90.0 ,arr.ind = TRUE))
AllUMAlpha<-as.vector(which(RNAseq_Ave$MAlphapercents >= 90.0 ,arr.ind = TRUE))
AllURSVAlpha<-as.vector(which(RNAseq_Ave$RSVAlphapercents >= 90.0 ,arr.ind = TRUE))

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

