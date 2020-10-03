{
  {
    
    sample<-"vivo_rep1"
    load(sprintf("%s/midif_pol_%s.RData",sample,sample))
    midif_pol2_vivo_rep1<-midif_pol2_dat
    
    sample<-"vivo_rep2"
    load(sprintf("%s/midif_pol_%s.RData",sample,sample))
    midif_pol2_vivo_rep2<-midif_pol2_dat
    
    #merge
    #vivo
    midif_pol_vivo<-rbind(midif_pol2_vivo_rep1,midif_pol2_vivo_rep1)
    ####
    gtf<-read_tsv("~/disk3/Nascent_RNA/reference/yeast/Ensembl/Saccharomyces_cerevisiae.R64-1-1.95.gtf",
                  col_names = F,comment = "#")%>%filter(X3=="gene")%>%
      mutate(len=X5-X4+1)%>%dplyr::select(gene=X9,len)
    gtf$gene <- str_split_fixed(gtf$gene,";",n=Inf)[,1] %>%str_sub(start = 9)%>%gsub('\\"',"",.) 
  }
  
  #############
  mysample<-"vivo"
  
  gene_type<-read_tsv("gene_type.txt",col_names = T) ;
  colnames(gene_type)=c("gene","type")
  mRNA<-gene_type%>%filter(type=="protein_coding")%>%getElement("gene")
  
  allGene<-intersect((eval(parse(text = sprintf("midif_pol_%s",mysample))))$gene.x%>%unique,mRNA) 
  allGene<-(eval(parse(text = sprintf("midif_pol_%s",mysample))))$gene.x%>%unique

  pol2_density<-mclapply(allGene,function(thisGene){
    len<-gtf%>%filter(gene==thisGene)%>%getElement("len")
    if(length(len)==0){return(NULL)}
    eachGene<- eval(parse(text = sprintf("midif_pol_%s",mysample)))%>%filter(gene.x==thisGene)%>%
      filter(gene.x==gene.y)
    footden<-eachGene%>%group_by(pol_pos)%>%dplyr::summarise(reads=length(pol_pos))
    mydata.fame<-data.frame(base=1:len,reads=0)
    mydata.fame$reads<-footden$reads[match(mydata.fame$base,footden$pol_pos)]
    mydata.fame[is.na(mydata.fame)] <- 0
    normalized_value<- mean(mydata.fame$reads,na.rm = T)
    #mean_perBase<-sum(mydata.fame$reads)/len
    if(normalized_value==0){return(NULL)}
    if(normalized_value >0){
      mydata.fame$reads <-mydata.fame$reads/normalized_value
      return(mydata.fame)  
    }
  },mc.cores = 40)%>%rbind.fill()
  
  
  
  dfdat<-pol2_density%>%group_by(base)%>%
    dplyr::summarise(dfdensity=mean(reads,na.rm=T))
  
  
  suppl.fig2a<- ggplot(dfdat%>%filter(base %in% c(1:1000) ),aes(x=base,y=dfdensity))+geom_point()+
    mystyle.print()+
    xlab("Distance from TSS")+
    ylab("Mean converage of reverse\nreads(transcription site)")
  
}