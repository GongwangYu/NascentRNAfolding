library(plyr)
library(dplyr)
library(readr)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(scales)
library(parallel)
#install.packages("DescTools")
library(DescTools)
options(stringsAsFactors =F )
{
  
  
  sample<-"vivo_rep1"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_rep1<-midif_pol2_dat
  
  sample<-"vivo_rep2"
  load(sprintf("%s/midif_pol_%s.RData",sample,sample))
  midif_pol2_rep2<-midif_pol2_dat
  
  #merge
  #vivo
  midif_pol_vivo<-rbind(midif_pol2_rep1,midif_pol2_rep2)
  ####
  
  ####
  gtf<-read_tsv("~/disk3/Nascent_RNA/reference/yeast/Ensembl/Saccharomyces_cerevisiae.R64-1-1.95.gtf",
                col_names = F,comment = "#")%>%filter(X3=="gene")%>%
    mutate(len=X5-X4+1)%>%dplyr::select(gene=X9,len)
  gtf$gene <- str_split_fixed(gtf$gene,";",n=Inf)[,1] %>%str_sub(start = 9)%>%gsub('\\"',"",.) 
  
}

###
#fig2a
{

  thisGene<-"RDN5-1";gene_len<-121 
  RDN<-midif_pol_vivo%>%filter(gene.x==thisGene)%>%filter(gene.x==gene.y)%>%
    filter(pol_pos<=gene_len)%>%mutate(rt_pos=rt_pos-1)%>%filter(rt_pos>0)
  
  RDN<-RDN%>%dplyr::group_by(pol_pos)%>%dplyr::mutate(num_rt=length(unique(rt_pos)))%>%filter(num_rt>2)%>%
    ungroup(pol_pos)
  
  pol2end = sort(RDN$pol_pos%>%unique())
  
  dfdat<-lapply(pol2end, function(i){
    eachrow<-RDN%>%filter(pol_pos==i)%>% dplyr::select(base=rt_pos)%>%
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
  rownames(dfdat)<-pol2end ; colnames(dfdat)<-1:gene_len

  rem_row=12
  toPlot <- melt(dfdat[-(1:rem_row),]%>%as.matrix()) 
  colnames(toPlot) <- c("Transcript length (nt)","Nucleotide position (nt)","Density of RT stop");

fig2a<-  toPlot %>% 
    ggplot(aes(y=`Transcript length (nt)`,x=`Nucleotide position (nt)`,
               fill=`Density of RT stop`)) +
    geom_tile(color= NA) +
    scale_y_reverse(expand=c(0,0), 
                    breaks =seq(40,120,10) ) +
    scale_x_continuous(expand = c(0,0), breaks =c(1,seq(0,gene_len,10)[-1]) )+
    theme_classic()+
    scale_fill_gradient("Density of RT stop",
                        #colours = rev(brewer.pal(11,"Spectral")),na.value=NA,
                        na.value=NA, low = "grey", high = "blue",
                        #breaks=c(-3,-7,-11),labels=expression(10^-3,10^-7,10^-11),
                        guide=guide_colorbar(title.position ="top" ,direction="horizontal",title.hjust=0.5,
                                             barwidth =8,barheight = 1, ticks.colour="black",ticks.linewidth = 1 ))+
    theme(legend.position = "top")
    
##################

plot.dat<-data.frame(base=pol2end,react=dfdat$`42`)%>%na.omit()

fig2c<- ggplot(plot.dat,aes(x=base,y=react))+geom_line()+geom_point()+
  xlab("Transcript length (nt)")+
  ylab("Density of RT stop" )+
  mystyle.print()
}

#########################
#fig2d
{




dffoldingChange<-mclapply(midif_pol_vivo$gene.x%>%unique,function(thisGene){
  
  len<-gtf%>%filter(gene==thisGene)%>%getElement("len")
  
  dfraw<-midif_pol_vivo%>%filter(gene.x==thisGene)%>%filter(gene.x==gene.y)%>%
    filter(pol_pos<=len)%>%mutate(rt_pos=rt_pos-1)%>%filter(rt_pos>0)
  
  dfraw<-dfraw%>%dplyr::group_by(pol_pos)%>%
    dplyr::mutate(num_rt=length(unique(rt_pos)))%>%filter(num_rt>3)%>%
    ungroup(pol_pos)
  if(nrow(dfraw)==0){return(NULL)}
  
  pol2end = sort(dfraw$pol_pos%>%unique())
  
  dfdat<-mclapply(pol2end, function(i){
    eachrow<-dfraw%>%filter(pol_pos==i)%>% dplyr::select(base=rt_pos)%>%
      group_by(base)%>%dplyr::summarise(read=length(base))
    eachdat<-data.frame(base=1:i)
    #df.each.dat<-full_join(eachdat,eachrow)
    df.each.dat<-merge(eachdat,eachrow,all = T)
    df.each.dat[is.na(df.each.dat)]<-NaN
    df.each.dat$read<-DescTools::Winsorize(df.each.dat$read, probs = c(0.05, 0.95),na.rm = T)
    df.each.dat$read<-df.each.dat$read/max(df.each.dat$read,na.rm = T)
    mydata.fame<-data.frame(base=1:len)
    # dfdat<-full_join(mydata.fame,df.each.dat)%>%dplyr::select(read)%>%t()%>%as.data.frame()
    dfdat<-merge(mydata.fame,df.each.dat,all=T)%>%dplyr::select(read)%>%t()%>%as.data.frame()
    return(dfdat)
  },mc.cores = 20)%>%rbind.fill()
  
  rownames(dfdat)<-pol2end ;  colnames(dfdat)<-1:len
  #dfdat<-dfdat%>%filter(`1`!= "NaN")
  if(nrow(dfdat)<50){return(NULL)}
  
  base=1
  foldingChange<- mclapply(1:len,function(base){
    dat<-dfdat[,base]%>%as.data.frame()%>%filter(.!="NA")%>%getElement(".")
    if(length(dat)<50) {return(NULL)}  
    
    window=50;step=1
    numberWin<-((length(dat)-window)%/%step)+1
    
    df_pvalue<-lapply(0:(numberWin-1),function(i){
      eachWinIndex<-(1+step*i):(window+step*i)
      eachwin<-dat[(1+step*i):(window+step*i)]
      myleft=eachwin[1:(length(eachwin)/2)]
      myright=eachwin[(length(eachwin)/2+1):length(eachwin)]
      if(length(na.omit(myleft))==0){return(NULL)}
      if(length(na.omit(myright))==0){return(NULL)}
      
      mytest<- wilcox.test(myleft,myright)
      myratio<-mean(myleft,na.rm = T)/mean(myright,na.rm = T)
      return(data.frame(pvalue=mytest$p.value,myratio))
    })%>%rbind.fill()
    
    if(is.null(df_pvalue)){return(NULL)}
    
    if(min(df_pvalue$pvalue,na.rm=T)<0.05){
      df_min_p<-df_pvalue%>%filter(pvalue==min(df_pvalue$pvalue,na.rm=T)) %>%
        getElement("myratio")
      return(data.frame(base=base,change=T,dfdir=df_min_p,gene=thisGene,len))
    } else {
      return(data.frame(base=base,change=F,dfdir=0,gene=thisGene,len))
    }
    
  },mc.cores = 20)%>%rbind.fill()
  
  return(foldingChange)
},mc.cores = 40)%>%rbind.fill()


save(dffoldingChange,file="./RData.file/dffoldingChange_final.RData")


dfper<-dffoldingChange%>%group_by(gene)%>%dplyr::summarise(per=sum(change==TRUE)/length(change))
gene_type<-read_tsv("gene_type.txt",col_names = T) ;
colnames(gene_type)=c("gene","type")
dfper<-dfper%>%merge(gene_type)

dfper_plot<-data.frame(per=dfper%>%filter(type=="protein_coding")%>%getElement("per"),
                       type="mRNA")%>%
  rbind(data.frame(per=dfper%>%filter(type!="protein_coding")%>%
                     getElement("per"),type="ncRNA"))


fig2d<-ggplot(dfper_plot,aes(x=per,y= ..count../sum(..count..),fill=type))+
  geom_histogram(data = subset(dfper_plot,type=="mRNA"), alpha = 0.7)+
  geom_histogram(data = subset(dfper_plot,type=="ncRNA"), alpha = 0.7)+
  scale_color_manual(values = c("red","blue"))+
  scale_x_continuous(breaks = c(0,0.25,0.50,0.75,1.00),
                     labels = c("0%","25%","50%","75%","100%"))+
  xlab("Percentage of bases with structural transition ")+
  ylab("Fraction of genes")+
  mystyle.print()+
  theme(legend.position = c(0.7,0.7),
        legend.title = element_blank())

ggarrange(fig2a,"", fig2c,fig2d,labels = c("a","b","c","d"),nrow = 2,ncol = 2,
          heights = c(2:1))
}
