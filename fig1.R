#####
#fig1b

{
  sample<-"vivo_rep1"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_rep1<-midif_pol2_dat
  
  sample<-"vivo_rep2"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_rep2<-midif_pol2_dat
  
  rep_1<-midif_pol2_rep1 %>%group_by(gene.x)%>%
    dplyr::summarise(rep1=length(gene.x)*10^6/nrow(midif_pol2_rep1))
  
  rep_2<-midif_pol2_rep2 %>%group_by(gene.x)%>%
    dplyr::summarise(rep2=length(gene.x)*10^6/nrow(midif_pol2_rep2))
  
merge_rep1_rep2<- merge(rep_1,rep_2)
  
cor.test(merge_rep1_rep2$rep1,merge_rep1_rep2$rep2,method = "p")
  #rho=0.99
  
fig1b<-  ggplot(merge_rep1_rep2,aes(x=rep_1,y=rep_2))+geom_point()+
  geom_abline(slope = 1,intercept = 0,lwd=1.5,col="red",lty=5)+
    scale_x_log10(limits=c(0.1,10000), breaks=c(1,100,10000),
                  labels= as.character(c(1,100,10000)))+
    scale_y_log10(limits=c(0.1,10000),breaks=c(1,100,10000),
                  labels= as.character(c(1,100,10000)))+
  xlab(label ="Reads per million (replicate 1)")+
    ylab("Reads per million(replicate 2)")+
    annotate("text",x=10,y=6000,label=expression("R = 0.99 P >" *10^-300))+
    mystyle.print()
 
} 

#######
#fig1c/d
{
gtf<-read_tsv("~/disk3/Nascent_RNA/reference/yeast/Ensembl/Saccharomyces_cerevisiae.R64-1-1.95.gtf",
              col_names = F,comment = "#")%>%dplyr::select(X3,X4,X5,X9)%>%filter(X3=="gene")
colnames(gtf)<-c("gene_type","start","end","gene")
gtf$gene <- str_split_fixed(gtf$gene,";",n=Inf)[,1] %>%str_sub(start = 9)%>%gsub('\\"',"",.) 


allGene<-intersect((midif_pol2_rep1$gene.x%>%unique()),(midif_pol2_rep2$gene.x%>%unique()))

modif_eachbase_rep1<-mclapply(allGene,function(thisGene){
  dfgtf<-gtf%>%filter(gene==thisGene) ; if(nrow(gtf)==0){return(NULL)}
  len=dfgtf$end-dfgtf$start+1
  eachGene<-midif_pol2_rep1%>%filter(gene.x==thisGene)%>%
    filter(gene.x==gene.y)%>%filter(rt_pos<=len &pol_pos<=len)
  df.rt<-eachGene%>%group_by(base=rt_pos)%>%dplyr::summarise(rt_num=length(rt_pos))
  df.pol2<-eachGene%>%group_by(base=pol_pos)%>%dplyr::summarise(pol2_num=length(pol_pos))
  mydata.fame<-data.frame(base=1:len,gene=thisGene)
  dfdat<-merge(mydata.fame,df.rt,all = T)%>%merge(df.pol2,all=T)
  return(dfdat)
},mc.cores = 40)%>%rbind.fill()

modif_eachbase_rep2<-mclapply(allGene,function(thisGene){
  dfgtf<-gtf%>%filter(gene==thisGene) ; if(nrow(gtf)==0){return(NULL)}
  len=dfgtf$end-dfgtf$start+1
  eachGene<-midif_pol2_rep2%>%filter(gene.x==thisGene)%>%
    filter(gene.x==gene.y)%>%filter(rt_pos<=len &pol_pos<=len)
  df.rt<-eachGene%>%group_by(base=rt_pos)%>%dplyr::summarise(rt_num=length(rt_pos))
  df.pol2<-eachGene%>%group_by(base=pol_pos)%>%dplyr::summarise(pol2_num=length(pol_pos))
  mydata.fame<-data.frame(base=1:len,gene=thisGene)
  dfdat<-merge(mydata.fame,df.rt,all = T)%>%merge(df.pol2,all=T)
  return(dfdat)
},mc.cores = 40)%>%rbind.fill()


allGene1<-intersect(modif_eachbase_rep1$gene%>%unique,modif_eachbase_rep2$gene%>%unique)

dfcor_rep_rt<-mclapply(allGene1,function(thisGene){
  myrep1<-modif_eachbase_rep1%>%filter(gene==thisGene)
  myrep2<-modif_eachbase_rep2%>%filter(gene==thisGene)
  dfmerge<-merge(myrep1,myrep2,by="base")
  dfmerge[is.na(dfmerge)]<-0
  if(nrow(filter(dfmerge,rt_num.x !=0 & rt_num.y !=0))<=10){return(NULL)}
  rt_average<-min(mean(dfmerge$rt_num.x,na.rm =T),mean(dfmerge$rt_num.y,na.rm = T))
  rt_cor<-cor.test(~rt_num.x+rt_num.y,data=dfmerge,method = "p")
  data.frame(gene=thisGene,
             rt_average,rt_cor=rt_cor$estimate)
},mc.cores = 40)%>%rbind.fill()

dfcor_rep_pol2<-mclapply(allGene1,function(thisGene){
  myrep1<-modif_eachbase_rep1%>%filter(gene==thisGene)
  myrep2<-modif_eachbase_rep2%>%filter(gene==thisGene)
  dfmerge<-merge(myrep1,myrep2,by="base")
  dfmerge[is.na(dfmerge)]<-0
  if(nrow(filter(dfmerge,pol2_num.x !=0 & pol2_num.y !=0))<=10){return(NULL)}
  pol2_average<-min(mean(dfmerge$pol2_num.x,na.rm =T),mean(dfmerge$pol2_num.y,na.rm = T))
  pol2_cor<-cor.test(~pol2_num.x+pol2_num.y,data=dfmerge,method = "p")
  data.frame(gene=thisGene,
             pol2_average,pol2_cor=pol2_cor$estimate
  )
  
},mc.cores = 40)%>%rbind.fill()

dfcor_rep_pol2<-dfcor_rep_pol2%>%filter(pol2_average>=1)

fig1c<-ggplot(dfcor_rep_rt,aes(x=rt_average,y=rt_cor))+geom_point()+
  geom_point(data=dfcor_rep_rt%>%filter(gene=="RDN5-1"),
             aes(x=rt_average,y=rt_cor),colour="red")+
  scale_x_log10(limits=c(1,1000),breaks = c(1,10,100,1000),labels =as.character(c(1,10,100,1000)))+
  scale_color_manual(values = c("black","red"))+
  xlab("Average number of reads of the gene")+
  ylab("Pearson correlation of coverage of\nforward reads between replicates")+
  theme(legend.position = "")+
  mystyle.print()+
  ylim(c(0.4,1))



fig1d<-ggplot(dfcor_rep_pol2,aes(x=pol2_average,y=pol2_cor))+geom_point()+
  geom_point(data=dfcor_rep_pol2%>%filter(gene=="RDN5-1"),
             aes(x=pol2_average,y=pol2_cor),colour="red")+
  scale_x_log10(limits=c(1,1000),breaks = c(1,10,100,1000),labels =as.character(c(1,10,100,1000)))+
  scale_color_manual(values = c("red","black"))+
  xlab("Average number of reads of the gene")+
  ylab("Pearson correlation of coverage of\nreverse reads between replicates")+
  theme(legend.position = "")+
  mystyle.print()+
  ylim(c(0.6,1))

##################
rep1<-modif_eachbase_rep1%>%filter(gene== "RDN5-1")
rep1[is.na(rep1)]<-0 ; rep1$type<-factor("rep1")

rep2<-modif_eachbase_rep2%>%filter(gene== "RDN5-1")
rep2[is.na(rep2)]<-0 ; rep2$type<-factor("rep2")

dfrep1_rep2<-rbind(rep1,rep2)


fig1c.inset<-ggplot(dfrep1_rep2,aes(x=base,y=rt_num,fill=type))+
  geom_bar(stat = "identity",binwidth = 0.01)+
  facet_grid(type~.) +
  ylab("Coverage of forward reads\n(RT stop)")+
  scale_x_continuous("Nucleotide positions(nt)", breaks = c(1,40,80,120))+
  mystyle.print()+
  theme(strip.text = element_blank(),
        legend.position = "none")
fig1c.inset

fig1d.inset<-ggplot(dfrep1_rep2,aes(x=base,y=pol2_num,fill=type))+
  geom_histogram(stat = "identity")+
  facet_grid(type~.) +
  ylab("Coverage of reverse reads\n( transcription site )")+
  scale_x_continuous("Nucleotide positions(nt)",breaks = c(1,40,80,120))+
  mystyle.print()+
  theme(strip.text = element_blank(),
        legend.position = "none")
fig1d.inset
}

#########################
#fig1e
{
  
  
  #merge
  #vivo
  midif_pol_vivo<-rbind(midif_pol2_rep1,midif_pol2_rep2)
  ####
    
  
  thisGene<-"tE(UUC)Q";gene_len<-72
  
  tRNA<-midif_pol_vivo%>%filter(gene.x==thisGene)%>%filter(gene.x==gene.y)%>%
    filter(pol_pos<=gene_len)%>%mutate(rt_pos=rt_pos-1)%>%filter(rt_pos>0)
  
  tRNA<-tRNA%>%dplyr::group_by(pol_pos)%>%dplyr::mutate(num_rt=length(unique(rt_pos)))%>%filter(num_rt>2)%>%
    ungroup(pol_pos)
  
  pol2end = sort(tRNA$pol_pos%>%unique())
  
  dfdat<-lapply(pol2end, function(i){
    eachrow<-tRNA%>%filter(pol_pos==i)%>% dplyr::select(base=rt_pos)%>%
      group_by(base)%>%dplyr::summarise(read=length(base))
    eachdat<-data.frame(base=1:i)
    df.each.dat<-full_join(eachdat,eachrow)
    df.each.dat[is.na(df.each.dat)]<-0
    df.each.dat$read<-DescTools::Winsorize(df.each.dat$read, probs = c(0.05, 0.95))
    df.each.dat$read<-df.each.dat$read/max(df.each.dat$read)
    mydata.fame<-data.frame(base=1:gene_len)
    dfdat<-full_join(mydata.fame,df.each.dat)%>%dplyr::select(read)%>%t()%>%as.data.frame()
    return(dfdat)
  })%>%rbind.fill()
  
  rownames(dfdat)<-pol2end
  colnames(dfdat)<-1:gene_len
  
  rem_row=9
  
  dfdat<- dfdat[-(1:rem_row),]
  dfdat<-dfdat[1:11,]
  toPlot <- melt(dfdat%>%as.matrix()) ;
  
  
  colnames(toPlot) <- c("Transcript length (nt)","Nucleotide position (nt)","Density of RT stop");
fig1e<-  toPlot %>% 
    ggplot(aes(y=`Transcript length (nt)`,x=`Nucleotide position (nt)`,
               fill=`Density of RT stop`)) +
    geom_tile(color= NA) +
    #scale_fill_gradient(na.value=NA, low = "grey", high = "blue") + 
    scale_y_reverse(expand=c(0,0),  breaks =seq(min(pol2end[-c(1:rem_row)]),65,1) ) +
    scale_x_continuous(expand = c(0,0), breaks =c(1,seq(0,gene_len,10)[-1]) )+
    theme_classic()+
    scale_fill_gradient("Density of RT stop",
                        #colours = rev(brewer.pal(11,"Spectral")),na.value=NA,
                        na.value=NA, low = "grey", high = "blue",
                        #breaks=c(-3,-7,-11),labels=expression(10^-3,10^-7,10^-11),
                        guide=guide_colorbar(title.position = "bottom",title.hjust=0.5,
                                            ticks.colour="black",ticks.linewidth = 1 ))+
    theme(legend.position = "bottom")
  fig1e
}



p1=ggarrange(fig1b,fig1c,fig1c.inset,nrow = 1,ncol = 3,labels = c("b","c"))

p2=ggarrange(fig1e,ggarrange(fig1d,fig1d.inset,"","",ncol = 2,nrow = 2),
             ncol = 2,labels = c("e","d"),widths = c(1,2))


ggarrange(p1,p2,nrow = 2,heights = c(1,2))
