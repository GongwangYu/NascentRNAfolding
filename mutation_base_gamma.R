gtf<-read_tsv("/mnt/data/home/yugongwang/reference/yeast/Ensembl/Saccharomyces_cerevisiae.R64-1-1.95.gff3",
              col_names = F,comment = "#")%>%dplyr::select(X1,X3,X4,X5,X7,X9)%>%
              filter(X3 %in% c("gene","ncRNA_gene"))
colnames(gtf)<-c("chr","gene_type","start","end","strand","gene")

gtf$gene <- str_split_fixed(gtf$gene,";",n=Inf)[,1] %>% str_sub(start = 9) 

mutation_perSite_raw<-read_tsv("mutationRate_perSite_gamma",col_names = T) 

heterozygous<-read.table("heterozygous_revised.tsv",sep = " ",header = T,stringsAsFactors = F)
colnames(heterozygous)[1]<-"chr_mut"

heterozygous<- heterozygous%>%mutate(chr=as.numeric(gsub("chr","",str_split_fixed(heterozygous$chr_mut,"_",2)[,1])),
                                     pos=as.numeric(str_split_fixed(heterozygous$chr_mut,"_",2)[,2]))

heterozygous$chr%>%unique()


xinaLuoma<-data.frame(xina=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"),
                      luoma=c(1:16))

heterozygous$chr<-xinaLuoma$xina[match(heterozygous$chr,xinaLuoma$luoma)]


mutation_perSite<- mclapply(1:421773,function(i){
  eachrow<-mutation_perSite_raw[i,] 
  chr<-heterozygous[i,193]
  pos<-heterozygous[i,194]
  numb<- strsplit(eachrow$`site  mp  m  k`, " ", fixed=TRUE) %>%unlist()
  return(data.frame(chr=chr,pos=pos,
                    mut=numb[length(numb)]%>%as.numeric()))
  
},mc.cores = 40)%>%rbind.fill()



thisGene<-gtf$gene[2]

dfmutation_eachgene<- mclapply(gtf$gene,function(thisGene){
  gene<-gtf%>%filter(gene==thisGene)
  len<-gene$end-gene$start+1
  if(gene$strand=="+"){
    eachGene_muat<-mutation_perSite%>%filter(chr==gene$chr) %>%filter(pos %in% c(gene$start:gene$end)) %>%
      mutate(base=pos-gene$start+1)
    if(nrow(eachGene_muat)==0){mydata.fame<-data.frame(base=1:len,mutRate=0,gene=thisGene)}else{
      mydata.fame<-data.frame(base=1:len,mutRate=0,gene=thisGene)
      mydata.fame$mutRate<-eachGene_muat$mut[match(mydata.fame$base,eachGene_muat$base)]
      mydata.fame[is.na(mydata.fame)] <- 0
       
    }
    return(mydata.fame)
  }
  
  if(gene$strand=="-"){
    eachGene_muat<-mutation_perSite%>%filter(chr==gene$chr) %>%filter(pos %in% c(gene$start:gene$end)) %>%
      mutate(base=gene$end -pos  +1)
    if(nrow(eachGene_muat)==0){mydata.fame<-data.frame(base=1:len,mutRate=0,gene=thisGene)}else{
      mydata.fame<-data.frame(base=1:len,mutRate=0,gene=thisGene)
      mydata.fame$mutRate<-eachGene_muat$mut[match(mydata.fame$base,eachGene_muat$base)]
      mydata.fame[is.na(mydata.fame)] <- 0
      
    }
    return(mydata.fame)
  }
  
  
},mc.cores = 40)

names(dfmutation_eachgene)<-gtf$gene
save(dfmutation_eachgene,file="dfmutation_eachgene.RData")
