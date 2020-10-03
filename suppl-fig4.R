load("RData.file/modif_react_vivo_4.RData")
load("DMS_seq/mRNA_foldingStrength.RData")
load("mystyle.print.RData")
#表达量
dfExpr <- read_tsv("Expression_Yeast_Snyder_2008.txt",col_names=TRUE);
colnames(dfExpr) <- c("gene","expr");
#保守性
load("/mnt/data/home/phil/acting/morphVari/rawdata/Yang2015/01.yeast.RData")
conservation_yeast<-  data.frame(gene=toPlot.yeast$orf,
                                 conser=log(1/toPlot.yeast$dNdS))


###########
#suppl.fig4a
{
allGene<-union(names(modif_react),dfExpr$gene)%>%
  union(conservation_yeast$gene)%>%
  union(mRNA_foldingStrength$gene)

#thisGene<-allGene[100]
dfreact_adaptive<- mclapply(allGene,function(thisGene){
  react<- modif_react[[thisGene]]%>%getElement("RSR") %>%mean()
  react<-(-react)
  #react<-(1/react)
  react_gini<-modif_react[[thisGene]] %>%getElement("RSR") %>%Gini()
  
  mature_rna_react<-mRNA_foldingStrength%>%filter(gene==thisGene) %>%getElement("dfreact")
  mature_rna_react<-(-mature_rna_react)
  if(length(mature_rna_react)==0){mature_rna_react<-NaN}
  
  mature_rna_react_gini<-mRNA_foldingStrength%>%filter(gene==thisGene) %>%getElement("dfgini")
  if(length(mature_rna_react_gini)==0){mature_rna_react_gini<-NaN}
  
  expr<-dfExpr%>%filter(gene==thisGene) %>%getElement("expr")
  if(length(expr)==0){expr<-NaN}
  
  conser<-conservation_yeast%>%filter(gene==thisGene) %>%getElement("conser")
  if(length(conser)==0){conser<-NaN}
  
  return(data.frame(gene=thisGene,react=react,react_gini=react_gini,
                    mature_rna_react,mature_rna_react_gini,
                    expr=expr,conser=conser))
},mc.cores = 40)%>%rbind.fill()

#nascnet react and expression /control muture rna react
dfdat<-dfreact_adaptive %>%dplyr::select(react,expr,mature_rna_react)%>%na.omit()
pcor.test(dfdat$react,dfdat$expr,dfdat$mature_rna_react,method = "s")
#n=1360,rho=0.3168, p <10^-32

dfdat<-dfreact_adaptive %>%dplyr::select(react_gini,expr,mature_rna_react)%>%na.omit()
pcor.test(dfdat$react_gini,dfdat$expr,dfdat$mature_rna_react,method = "s")
#n=1360,rho=0.5616,p <10^-113

##nascnet react and conservation /control muture rna react
dfdat<-dfreact_adaptive %>%dplyr::select(react,conser,mature_rna_react)%>%na.omit()
pcor.test(dfdat$react,dfdat$conser,dfdat$mature_rna_react,method = "s")
#n=1285,rho=0.1795,p<10^-10

dfdat<-dfreact_adaptive %>%dplyr::select(react_gini,conser,mature_rna_react)%>%na.omit()
pcor.test(dfdat$react_gini,dfdat$conser,dfdat$mature_rna_react,method = "s")
#n=1285,rho=0.3784,p< 10^-44

dfpcor_plot<-data.frame(rho=c(0.317,0.562,0.180,0.378),
                        type=(c("Nascent RNA folding strength(Negative of \nsingle-stranded score)~expression level",
                                "Nascent RNA folding strength(Gini index of \nsingle-stranded score)~expression level",
                                "Nascent RNA folding strength(Negative of \nsingle-stranded score)~evolutionary conservation",
                                "Nascent RNA folding strength(Gini index of \nsingle-stranded score)~evolutionary conservation")),
                        id=c(4:1))

suppl.fig4a<-ggplot(dfpcor_plot,aes(x=id,y=rho))+geom_bar(stat="identity")+
  coord_flip()+
  scale_x_continuous(breaks =dfpcor_plot$id,labels=dfpcor_plot$type)+
  ylab("Partial correlation controlling mature RNA structure")+xlab("")+
  mystyle.print()+
  theme(axis.text.y =element_text(size = unit(12,"bigpts")))+
  annotate("text",x=4,y=0.7,label=expression(italic(P) * " < " * 10^{-32}))+
  annotate("text",x=3,y=0.7,label=expression(italic(P) * "< " * 10^{-113}))+
  annotate("text",x=2,y=0.7,label=expression(italic(P) * " < " * 10^{-10}))+
  annotate("text",x=1,y=0.7,label=expression(italic(P) * " < " * 10^{-44}))+
  ylim(c(0,0.8))

suppl.fig4a
}
################
#suppl.fig4b
{
allGene<-union(names(modif_react),dfExpr$gene)%>%
  union(conservation_yeast$gene)

dfreact_adaptive<- mclapply(allGene,function(thisGene){
  react<- modif_react[[thisGene]]%>%getElement("RSR") %>%mean()
  react<-(-react)
  react_gini<-modif_react[[thisGene]] %>%getElement("RSR") %>%Gini()
  expr<-dfExpr%>%filter(gene==thisGene) %>%getElement("expr")
  if(length(expr)==0){expr<-NaN}
  conser<-conservation_yeast%>%filter(gene==thisGene) %>%getElement("conser")
  if(length(conser)==0){conser<-NaN}
  return(data.frame(gene=thisGene,react=react,react_gini=react_gini,
                    expr=expr,conser=conser))
},mc.cores = 40)%>%rbind.fill()



myrange<- quantile(mRNA_foldingStrength$dfreact,p=c(0.4,0.6))
contrled.gene<- mRNA_foldingStrength%>%filter(dfreact >=myrange[1]& dfreact <=myrange[2])%>%
  getElement("gene")%>%as.character()

cor.test(~react+expr,dfreact_adaptive%>%filter(gene %in% contrled.gene),method="s")
#N=1014 rho=0.246,p <10^-4 ,react_control

cor.test(~react_gini+expr,dfreact_adaptive%>%filter(gene %in% contrled.gene),method="s")
#N=1014,rho=0.637,p<10^-31,react_control

cor.test(~react+conser,dfreact_adaptive%>%filter(gene %in% contrled.gene),method="s")
#N=1014,rho=0.250,p < 10^-4 ,react_control

cor.test(~react_gini+conser,dfreact_adaptive%>%filter(gene %in% contrled.gene),method="s")
#N=1014,rho=0.350,p<10^-7,react_control

dfpcor_plot<-data.frame(rho=c(0.246,0.637,0.250,0.350),
                        type=(c("Nascent RNA folding strength(Negative of \nsingle-stranded score)~expression level",
                                "Nascent RNA folding strength(Gini index of \nsingle-stranded score)~expression level",
                                "Nascent RNA folding strength(Negative of \nsingle-stranded score)~evolutionary conservation",
                                "Nascent RNA folding strength(Gini index of \nsingle-stranded score)~evolutionary conservation")),
                        id=c(4:1))

suppl.fig4b<-ggplot(dfpcor_plot,aes(x=id,y=rho))+geom_bar(stat="identity")+
  coord_flip()+
  scale_x_continuous(breaks =dfpcor_plot$id,labels=dfpcor_plot$type)+
  ylab("Spearman's rank correlation among genes\n with similar mature RNA structure")+xlab("")+
  mystyle.print()+
  theme(axis.text.y =element_text(size = unit(12,"bigpts")))+
  annotate("text",x=4,y=0.8,label=expression(italic(P)* " < " * 10^{-4}))+
  annotate("text",x=3,y=0.8,label=expression(italic(P)* " < " * 10^-31))+
  annotate("text",x=2,y=0.8,label=expression(italic(P)*  " < " * 10^{-4}))+
  annotate("text",x=1,y=0.8,label=expression(italic(P)*  " < " * 10^{-7}))+
  ylim(c(0,0.9))

suppl.fig4b
}

ggarrange(suppl.fig4a,suppl.fig4b,nrow = 2,labels = c("a","b"))

