#TF Binding Site

library(enrichR)
dbs<-c("TRANSFAC_and_JASPAR_PWMs","ENCODE_TF_ChIP-seq_2015","Transcription_Factor_PPIs")
TBX21<-enrichr(genes = c("TBX21"),databases = dbs)
TBX21_JASPR<-TBX21$TRANSFAC_and_JASPAR_PWMs
#None of obvious interest
TBX21_ChIP<-TBX21$`ENCODE_TF_ChIP-seq_2015`
#None of obvious interest
TBX21_PPI<-TBX21$Transcription_Factor_PPIs
#GATA3

TBX21_AIM2<-enrichr(genes = c("TBX21","AIM2"),databases = dbs)
TBX21_AIM2_JASPR<-TBX21_AIM2$TRANSFAC_and_JASPAR_PWMs
#NFkB, STAT3, GATA1
TBX21_AIM2_ChIP<-TBX21_AIM2$`ENCODE_TF_ChIP-seq_2015`
#STAT1, STAT2, CTCF, IRF4, GATA2, GATA3
TBX21_AIM2_PPI<-TBX21_AIM2$Transcription_Factor_PPIs
#None of obvious interest

TBX21_ITK<-enrichr(genes = c("TBX21","ITK"),databases = dbs)
TBX21_ITK_JASPR<-TBX21_ITK$TRANSFAC_and_JASPAR_PWMs
TBX21_ITK_ChIP<-TBX21_ITK$`ENCODE_TF_ChIP-seq_2015`
TBX21_ITK_PPI<-TBX21_ITK$Transcription_Factor_PPIs
#None of obvious interest

library(ggplot2)
TargetGene_NOS2<-read.delim("NOS2.1.tsv")
write.csv(TargetGene_NOS2,"Complete Target Gene NOS2.csv")
subset_TargetGene_NOS2<-TargetGene_NOS2[c(1:50),]
subset_TargetGene_NOS2_Plot<-ggplot(subset_TargetGene_NOS2,aes(x=subset_TargetGene_NOS2$Target_genes,y=subset_TargetGene_NOS2$NOS2.Average, color=subset_TargetGene_NOS2$NOS2.Average))+geom_point()
subset_TargetGene_NOS2_Plot<-subset_TargetGene_NOS2_Plot+scale_color_gradient(name="NOS2 Average",low="blue",high="red")+theme_bw()
subset_TargetGene_NOS2_Plot<-subset_TargetGene_NOS2_Plot+labs(y="NOS2 Average",x="Target Gene")+coord_flip()
subset_TargetGene_NOS2_Plot<-subset_TargetGene_NOS2_Plot

TargetGene_TBX21<-read.delim("TBX21.1.tsv")
write.csv(TargetGene_TBX21,"Target Gene TBX21.csv")
subset_TargetGene_TBX21<-TargetGene_TBX21[c(1:50),]
subset_TargetGene_TBX21_Plot<-ggplot(subset_TargetGene_TBX21,aes(x=subset_TargetGene_TBX21$Target_genes,y=subset_TargetGene_TBX21$TBX21.Average, color=subset_TargetGene_TBX21$TBX21.Average))+geom_point()
subset_TargetGene_TBX21_Plot<-subset_TargetGene_TBX21_Plot+scale_color_gradient(name="TBX21 Average",low="blue",high="red")+theme_bw()
subset_TargetGene_TBX21_Plot<-subset_TargetGene_TBX21_Plot+labs(y="TBX21 Average",x="Target Gene")+coord_flip()
subset_TargetGene_TBX21_Plot<-subset_TargetGene_TBX21_Plot
#SLAMF8


Enrich_All_90<-read.delim("Enrichment_All_90.tsv",header=FALSE)
colnames(Enrich_All_90)<-c("ID","Antigen Class","Antigen","Cell Class","Cell","Num of peaks","Overlaps/ All 90 Perc Unique","Overlaps/ Control","Log P-value","Log Q-value","Fold enrichment","FE>1 ?")
Enrich_All_90_A549<-Enrich_All_90[Enrich_All_90$Cell == "A549",]

#No TFs overlap all genes from 90 perc unique

Enrich_CXCL11_CCL8<-read.delim("Enrichment CXCL11_CCL8.tsv",header=FALSE)
colnames(Enrich_CXCL11_CCL8)<-c("ID","Antigen Class","Antigen","Cell Class","Cell","Num of peaks","Overlaps/ CXCL11_CCL8","Overlaps/ Control","Log P-value","Log Q-value","Fold enrichment","FE>1 ?")
Enrich_CXCL11_CCL8_A549<-Enrich_CXCL11_CCL8[Enrich_CXCL11_CCL8$Cell == "A549",]
#No TFs overlap all genes from 90 perc unique

Enrich_TBX21<-read.delim("Enrichment_TBX21.tsv",header=FALSE)
colnames(Enrich_TBX21)<-c("ID","Antigen Class","Antigen","Cell Class","Cell","Num of peaks","Overlaps/ TBX21","Overlaps/ Control","Log P-value","Log Q-value","Fold enrichment","FE>1 ?")
Enrich_TBX21_A549<-Enrich_TBX21[Enrich_TBX21$Cell == "A549",]

