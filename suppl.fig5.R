load("RData.file/modif_react_vivo_4.RData")
load("mutation/dfmutation_eachgene.RData")
load("DMS_seq//mRNA_foldingStrength.RData")
load("R-loop/gene.r.score.RData")
load("mystyle.print.RData")
options(stringsAsFactors = F)


allGene<-union(names(modif_react),names(gene.r.score))%>%
  union(names(dfmutation_eachgene))%>%union(mRNA_foldingStrength$gene)

react_r.score_mut<-mclapply(allGene,function(thisGene){
  react<- modif_react[[thisGene]]%>%getElement("RSR") %>%mean()
  react<-(-react)
  
  react_gini<-modif_react[[thisGene]] %>%getElement("RSR") %>%Gini()
  
  mature_react<-mRNA_foldingStrength%>%filter(gene==thisGene) %>%getElement("dfreact")
  mature_react<-(-mature_react)
  if(length(mature_react)==0){mature_react<-NaN}
  
  # mature_react_gini<-mRNA_foldingStrength%>%filter(gene==thisGene) %>%getElement("dfgini")
  # if(length(mature_react_gini)==0){mature_react_gini<-NaN}
  
  mut<-dfmutation_eachgene[[thisGene]] %>%getElement("mutRate")%>%mean()
  
  if(is.null(gene.r.score[[thisGene]])){r.score<-NaN}else
  {r.score<-gene.r.score[[thisGene]] %>%getElement("r.score")}
  
  data.frame(gene=thisGene,mut,r.score,react, react_gini,mature_react)
},mc.cores = 40)%>%rbind.fill() %>%na.omit()


############suppl.fig5a 
dfdat<-react_r.score_mut 

#nascnet react and r-loop score /control muture rna react
pcor.test(dfdat$react,dfdat$r.score,dfdat$mature_react,method = "s")
#n=1351,rho= -0.022,p=0.4294896

pcor.test(dfdat$react_gini,dfdat$r.score,dfdat$mature_react,method = "s")
#n=1351,rho=-0.442,p<10^-65

##nascnet react and mutation rate /control muture rna react
pcor.test(dfdat$react,dfdat$mut,dfdat$mature_react,method = "s")
#n=1351,rho=-0.125,p<10^-5

pcor.test(dfdat$react_gini,dfdat$mut,dfdat$mature_react,method = "s")
#n=1351,rho=-0.219,p<10^-15



#################################
#suppl.fig5b
myrange<- quantile(react_r.score_mut$mature_react,p=c(0.4,0.6))

dfdat<-react_r.score_mut %>%filter(mature_react >=myrange[1] & mature_react <=myrange[2])

cor.test(~react+r.score,data = dfdat,method="s")
#n=271,rho=-0.057,p=0.34
cor.test(~react_gini+r.score,data = dfdat,method="s")
#n=271,rho=-0.411,p <10^-11
cor.test(~react+mut,data = dfdat,method="s")
#n=271,rho=-0.198,p<0.005
cor.test(~react_gini+mut,data = dfdat,method="s")
#n=271,rho=-0.244,p<10^-4


suppl.fig5a.dat<-data.frame(rho=c(-0.022,-0.442,-0.125,-0.219),
                            type=(c("Nascent RNA folding strength(Negative\n of single-stranded score)~R-loop score",
                                    "Nascent RNA folding strength(Gini index\n of single-stranded score)~R-loop score",
                                    "Nascent RNA folding strength(Negative\n of single-stranded score)~mutation rate",
                                    "Nascent RNA folding strength(Gini index\n of single-stranded score)~mutation rate")),
                            id=c(4:1))


suppl.fig5a<-  ggplot(suppl.fig5a.dat,aes(x=id,y=rho))+geom_bar(stat="identity")+
  coord_flip()+
  scale_x_continuous(breaks =suppl.fig5a.dat$id,labels=suppl.fig5a.dat$type)+
  ylab("Partial correlation controlling mature RNA structure")+xlab("")+
  mystyle.print()+
  ylim(c(-0.6,0))+
  theme(axis.text.y =element_text(size = unit(12,"bigpts")))+
  annotate("text",x=4,y=-0.5,label=expression(italic(P) * " = 0.43"))+
  annotate("text",x=3,y=-0.5,label=expression(italic(P) *" < " * 10^{-65}))+
  annotate("text",x=2,y=-0.5,label=expression(italic(P) * " < " *10^-5))+
  annotate("text",x=1,y=-0.5,label=expression(italic(P) * " < " * 10^{-15}))



suppl.fig5b.dat<-data.frame(rho=c(-0.057,-0.411,-0.198,-0.244),
                            type=(c("Nascent RNA folding strength(Negative\n of single-stranded score)~R-loop score",
                                    "Nascent RNA folding strength(Gini index\n of single-stranded score)~R-loop score",
                                    "Nascent RNA folding strength(Negative\n of single-stranded score)~mutation rate",
                                    "Nascent RNA folding strength(Gini index\n of single-stranded score)~mutation rate")),
                            id=c(4:1))


suppl.fig5b<-  ggplot(suppl.fig5b.dat,aes(x=id,y=rho))+geom_bar(stat="identity")+
  coord_flip()+
  scale_x_continuous(breaks =suppl.fig5b.dat$id,labels=suppl.fig5b.dat$type)+
  ylab("Spearman's rank correlation among genes\n with similar mature RNA structure")+xlab("")+
  mystyle.print()+
  ylim(c(-0.6,0))+
  theme(axis.text.y =element_text(size = unit(12,"bigpts")))+
  annotate("text",x=4,y=-0.5,label=expression(italic(P)* " = 0.34"))+
  annotate("text",x=3,y=-0.5,label=expression(italic(P)*" < " * 10^{-11}))+
  annotate("text",x=2,y=-0.5,label=expression(italic(P)*" < 0.005"))+
  annotate("text",x=1,y=-0.5,label=expression(italic(P)*" < " * 10^{-4}))





ggarrange(suppl.fig5a,suppl.fig5b,nrow = 2,labels = c("a","b"))


