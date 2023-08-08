
library(readr);
library(plyr);
library(dplyr);
library(parallel);
library(ggplot2);
library(Matrix);
library(ineq)
library(stringr)
############
{
  
  sample<-"vivo_rep1"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_vivo_rep1<-midif_pol2_dat
  
  sample<-"vivo_rep2"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_vivo_rep2<-midif_pol2_dat
  
  sample<-"dmso_rep1"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_dmso_rep1<-midif_pol2_dat
  
  sample<-"dmso_rep2"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_dmso_rep2<-midif_pol2_dat
  
  ####
  #merge
  #DMSO
  midif_pol_DMSO<-rbind(midif_pol2_dmso_rep1,midif_pol2_dmso_rep2)

  #vivo
  midif_pol_vivo<-rbind(midif_pol2_vivo_rep1,midif_pol2_vivo_rep1)
  
  ####
  gtf<-read_tsv("~/disk3/Nascent_RNA/reference/yeast/Ensembl/Saccharomyces_cerevisiae.R64-1-1.95.gtf",
                col_names = F,comment = "#")%>%filter(X3=="gene")%>%
    mutate(len=X5-X4+1)%>%dplyr::select(gene=X9,len)
  gtf$gene <- str_split_fixed(gtf$gene,";",n=Inf)[,1] %>%str_sub(start = 9)%>%gsub('\\"',"",.) 
}

#############
#Here, the threshold of RT stop/nt can be set to 1, 2, 4, in our analysis, we set this threshold to 1
coverage_threshold=1

for(mysample in c("vivo")){
  allGene<-intersect((eval(parse(text = sprintf("midif_pol_%s",mysample))))$gene.x%>%unique,
                     midif_pol_DMSO$gene.x%>%unique)
  
 
  Kdmso<-(nrow(midif_pol_DMSO)+nrow(eval(parse(text = sprintf("midif_pol_%s",mysample)))))/
    (2*nrow(midif_pol_DMSO))
  Ksample<-(nrow(midif_pol_DMSO)+nrow(eval(parse(text = sprintf("midif_pol_%s",mysample)))))/
    (2*nrow(eval(parse(text = sprintf("midif_pol_%s",mysample)))))

  modif_react<-mclapply(allGene,function(thisGene){
    len<-gtf%>%filter(gene==thisGene)%>%getElement("len")
    if(length(len)==0){return(NULL)}
    ###vivo
    eachGene<- eval(parse(text = sprintf("midif_pol_%s",mysample)))%>%filter(gene.x==thisGene)%>%
      filter(gene.x==gene.y)%>%group_by(base=rt_pos)%>%dplyr::summarise(reads=length(rt_pos))
    
    eachGene$base<-eachGene$base-1
    mydata.fame<-data.frame(base=0:len,reads=0)
    mydata.fame$reads<-eachGene$reads[match(mydata.fame$base,eachGene$base)]
    mydata.fame[is.na(mydata.fame)] <- 0
    rt.mean<-mean(mydata.fame$reads)
   
    if(rt.mean< coverage_threshold){return(NULL)}
    mydata.fame$react<-mydata.fame$reads*Ksample+1
    dfdat_sample<-mydata.fame%>%filter(base>0)
    
    ###DMSO
    eachGene<-midif_pol_DMSO%>%filter(gene.x==thisGene)%>%filter(gene.x==gene.y)%>%
      group_by(base=rt_pos)%>%dplyr::summarise(reads=length(rt_pos))
    
    eachGene$base<-eachGene$base-1
    mydata.fame<-data.frame(base=0:len,reads=0)
    mydata.fame$reads<-eachGene$reads[match(mydata.fame$base,eachGene$base)]
    mydata.fame[is.na(mydata.fame)] <- 0
    mydata.fame$react<-mydata.fame$reads*Kdmso+1
    dfdat_DMSO<-mydata.fame%>%filter(base>0)
    
    #raw structural reactivaty
    eachGene_RSR<- data.frame(
      base=1:len, 
      RSR=log2(dfdat_sample$react/dfdat_DMSO$react+1),

      gene=thisGene)
    
    eachGene_RSR$RSR<- ifelse(eachGene_RSR$RSR<7,eachGene_RSR$RSR,7)
   
    return(eachGene_RSR)
  },mc.cores = 40)
  
  names(modif_react)<-allGene
  save(modif_react,file = sprintf("RData.file/modif_react_%s_4.RData",mysample))
  
}

###################################################################end
