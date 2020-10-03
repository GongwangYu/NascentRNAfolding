
library(dplyr)
for(sample in c("vivo_rep1","vivo_rep2","dmso_rep1","dmso_rep2")){
  
  #forward read(RT stop)
  modif_dat<- read_tsv(sprintf("%s/%s_1.txt",sample,sample),col_names = F)%>%dplyr::select(X1,X3,X4)
  colnames(modif_dat)<-c("readID","gene","rt_pos")
  
  #reverse read(transcription site)
  Pol2_dat<-read_tsv(sprintf("%s/%s_2.txt",sample,sample),col_names = F)%>%dplyr::select(X1,X3,X4,X5)
  colnames(Pol2_dat) <-c("readID","gene","pos","cigar")
  
  Pol2_dat<-Pol2_dat%>%group_by(readID,gene)%>%
    dplyr::summarise(pol_pos=(str_extract_all(cigar,"\\d+",simplify=T) %>%as.numeric()%>%sum())+pos-1)
  
  
  midif_pol2_dat<-merge(modif_dat,Pol2_dat,by="readID") 
  
  save(midif_pol2_dat,file = sprintf("%s/midif_pol_%s.RData",sample,sample))
  
} 




