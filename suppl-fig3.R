

#suppl.fig3a
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

allGene<-(eval(parse(text = sprintf("midif_pol_%s",mysample))))$gene.x%>%unique

#only gene mean coverage >1
load("RData.file/modif_react_vivo_4.RData")
allGene<-names(modif_react)

pol2_footprint<-mclapply(allGene,function(thisGene){
  len<-gtf%>%filter(gene==thisGene)%>%getElement("len")
  if(length(len)==0){return(NULL)}
  if(is.null(modif_react[[thisGene]])){return(NULL)}
  
  eachGene<- eval(parse(text = sprintf("midif_pol_%s",mysample)))%>%filter(gene.x==thisGene)%>%
    filter(gene.x==gene.y)
  
  footprint<-eachGene%>%group_by(pol_pos)%>%filter(rt_pos==max(rt_pos))%>%
    dplyr::select(rt_pos,pol_pos)%>%unique()%>%dplyr::summarise(dfdistance=(pol_pos-rt_pos+1))
  return(footprint)
},mc.cores = 40)%>%rbind.fill()%>%filter(dfdistance >0)

save(pol2_footprint,file="RData.file/pol2_footprint.RData")


load("RData.file/pol2_footprint.RData")
load("mystyle.print.RData")
dfdat<-pol2_footprint%>%filter(dfdistance >14 &dfdistance<23)%>%group_by(dfdistance)%>%
  dplyr::summarise(dfcount=length(dfdistance))

suppl.fig3a<-ggplot(dfdat,aes(x=dfdistance,y=dfcount))+
  geom_point()+
  geom_line()+
  scale_x_continuous("Distance between transcription \nsite and proximal RT stop",
                     breaks =c(15:22) ,labels = c(15:22))+
  ylab("Number of read pairs ")+
  mystyle.print()

}

############
#suppl.fig3b
#{
thisGene<-"tE(UUC)Q";gene_len<-72#有用

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

toPlot <- melt(dfdat[-(1:rem_row),]%>%as.matrix()) ;


colnames(toPlot) <- c("Transcript length (nt)","Nucleotide position (nt)","single-stranded score");
toPlot %>% 
  ggplot(aes(y=`Transcript length (nt)`,x=`Nucleotide position (nt)`,
             fill=`single-stranded score`)) +
  geom_tile(color= NA) +
  #scale_fill_gradient(na.value=NA, low = "grey", high = "blue") + 
  scale_y_reverse(expand=c(0,0), breaks =seq(min(pol2end[-c(1:rem_row)]),max(pol2end),10) ) +
  scale_x_continuous(expand = c(0,0), breaks =c(1,seq(0,gene_len,10)[-1]) )+
  theme_classic()+
  scale_fill_gradient("single-stranded score",
                       #colours = rev(brewer.pal(11,"Spectral")),na.value=NA,
                       na.value=NA, low = "grey", high = "blue",
                       #breaks=c(-3,-7,-11),labels=expression(10^-3,10^-7,10^-11),
                       guide=guide_colorbar(title.position = "top",direction="horizontal",title.hjust=0.5,
                                            barwidth =10,ticks.colour="black",ticks.linewidth = 1 ))+
  theme(legend.position = "top")
