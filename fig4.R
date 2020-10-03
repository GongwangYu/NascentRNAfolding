###########
load("RData.file/modif_react_vivo_4.RData")

#表达量
dfExpr <- read_tsv("Expression_Yeast_Snyder_2008.txt",col_names=TRUE);
colnames(dfExpr) <- c("gene","expr");
#保守性
load("/mnt/data/home/phil/acting/morphVari/rawdata/Yang2015/01.yeast.RData")
conservation_yeast<-  data.frame(gene=toPlot.yeast$orf,
                                 conser=log(1/toPlot.yeast$dNdS))

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



cor.test(~react+expr,dfreact_adaptive,method="s")
#n=1404,rho=0.31,p<10e-32
cor.test(~react+conser,dfreact_adaptive,method="s")
#n=1311,rho=0.19,p<10e-11

cor.test(~react_gini+expr,dfreact_adaptive,method="s")
#n=1404,rho=0.56,p<10e-117
cor.test(~react_gini+conser,dfreact_adaptive,method="s")
#n=1311,rho=0.38,p<10e-46

#######

#画点图
{
  dfdat<-dfreact_adaptive%>%filter(expr!="NaN" & react!="NaN")
  fig3c<-ggplot(aes(x=expr,y=react),data = dfdat) +
    geom_point()+
    geom_smooth(method = 'lm',se=F,color='red',size=1)+
    scale_y_continuous("Nascent RNA folding strength\n(Negative of single-stranded score)") +
    scale_x_continuous("Expression level", breaks = c(0,5,10),labels=expression(2^0,2^5,2^10))+
    annotate("text",x=5,y=0,label=expression( italic(rho)* " = 0.31  " * italic(P) * " < " * 10^{-32}))+
    mystyle.print()
  
  dfdat<-dfreact_adaptive%>%filter(expr!="NaN" & react_gini!="NaN")
  fig3d<-ggplot(aes(x=expr,y=react_gini),data = dfdat) +
    geom_point()+
    geom_smooth(method = 'lm',se=F,color='red',size=1)+
    scale_y_continuous("Nascent RNA folding strength\n(Gini index of single-stranded score)") +
    scale_x_continuous("Expression level", breaks = c(0,5,10),labels=expression(2^0,2^5,2^10))+
    annotate("text",x=5,y=0.5,label=expression(italic(rho) * " = 0.56  " * italic(P)* " < " * 10^{-117}))+
    mystyle.print()
  
  dfdat<-dfreact_adaptive%>%filter(conser!="NaN" & react!="NaN")
  fig3e<-ggplot(aes(x=conser,y=react),data = dfdat) +
    geom_point()+
    geom_smooth(method = 'lm',se=F,color='red',size=1)+
    scale_y_continuous("Nascent RNA folding strength\n(Negative of single-stranded score)") +
    scale_x_continuous("Evolutionary conservation")+
    annotate("text",x=4,y=0,label=expression(italic(rho) *" = 0.19  " * italic(P) * " < " * 10^{-11}))+
    mystyle.print()
  
  dfdat<-dfreact_adaptive%>%filter(conser!="NaN" & react_gini!="NaN")
  fig3f<-ggplot(aes(x=conser,y=react_gini),data = dfdat) +
    geom_point()+
    geom_smooth(method = 'lm',se=F,color='red',size=1)+
    scale_y_continuous("Nascent RNA folding strength\n(Gini index of single-stranded score)") +
    scale_x_continuous("Evolutionary conservation")+
    annotate("text",x=5,y=0.5,label=expression(italic(rho)* " = 0.38  "* italic(P)* " < " * 10^{-46}))+
    mystyle.print()
  
}
