library(dplyr);
library(plyr);
library(ggplot2);
source("/home/phil/lib/R/phil/MyGgplot2Style.R");
load("mystyle.print.RData")

load("RData.file/modif_react_vivo_4.RData")
load("mutation/dfmutation_eachgene.RData")

allGene<-intersect(names(modif_react),names(dfmutation_eachgene))

mut_ratio<-mclapply(allGene,function(thisGene){
  eachreact<-modif_react[[thisGene]]
  if(is.null(eachreact)){return(NULL)}
  max_single<-order(eachreact$RSR,decreasing = T)[1:50]
  max_double<-order(eachreact$RSR,decreasing = F)[1:50]
  eachmut<-dfmutation_eachgene[[thisGene]]
  if(is.null(eachmut)){return(NULL)}
  if(length(eachmut$mutRate[eachmut$mutRate>0])/nrow(eachmut) < 0.05 ){return(NULL)}
  # if((nrow(eachmut%>%filter(base %in% max_single & mutRate >0)) <3) | 
  #    (nrow(eachmut%>%filter(base %in% max_double & mutRate >0)) <3) ){return(NULL)}
  
  #dfratio<-(mean(eachmut$mutRate[max_single])+0.001)/(mean(eachmut$mutRate[max_double])+0.001)
  dfratio<-mean(eachmut$mutRate[max_single])/(mean(eachmut$mutRate[max_double]))
  data.frame(gene=thisGene,dfratio)
},mc.cores = 20)%>%rbind.fill()%>%na.omit

mut_ratio$dfratio[mut_ratio$dfratio==Inf]<-50

dfExpr <- read_tsv("Expression_Yeast_Snyder_2008.txt",col_names=TRUE);
colnames(dfExpr) <- c("gene","expr");

dfmerge<-merge(mut_ratio,dfExpr)

cor.test(~dfratio+expr,data = dfmerge,method="s")
#N=112,rho=0.24,p<0.01

suppl.fig6a<-ggplot(dfmerge,aes(x=expr,y=dfratio))+geom_point()+
  scale_x_continuous("Expression level", breaks = c(0,3,6,9),labels=expression(2^0,2^3,2^6,2^9))+
  scale_y_continuous("Fold redution of mutation rate")+
  annotate("text",x=5,y=60,
           label=expression(paste(italic(rho) ," = 0.24   ", italic(P)," < ",0.01)))+
  geom_point(aes(x=x,y=y),fill="black",shape=15,data=data.frame(x=6.977280,y=7.26597848),size=3) +
  geom_point(aes(x=x,y=y),fill="black",shape=17,data=data.frame(x=5.149747,y=48.07214750),size=3) +
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  mystyle.print()

suppl.fig6a




###########################################
baselineMutRate <- 10 ^ seq(-10,-5,length.out=500); ## baseline mutation rate without the effect of nascent RNA folding
obsMutRate <- 10 ^ seq(-10,-5,length.out=500); ## observed mutation rate (have been reduced by nascent RNA folding)
foldReductionByNascentStr <- 10 ^ seq(0,2,length.out=500)[-1]; ## fold reduction of mutation rate due to nascent RNA folding
#fracDelMutation <- c(0.5,0.7,0.9,0.99,0.999,0.9999); ## fraction of deleterious mutations (this does not affect much)
fracDelMutation <- 0.7; ## fraction of deleterious mutations
affectRange <- 10; ## number of nucleotides affected by a single mutation


dfVal <- expand.grid(
  m = obsMutRate,
  r = foldReductionByNascentStr,
  d = fracDelMutation,
  l = affectRange) %>%
  mutate(s = m * (r - 1) * d * l);
dfLine.yeast <- dfVal %>%
  group_by(d,r) %>%
  arrange(s) %>%
  filter(s > 1e-7) %>%
  filter(row_number() == 1);
dfLine.yeast <- dfLine.yeast %>% rbind.fill(data.frame(m=1e-5,r=1,l=affectRange,d=fracDelMutation,s=0))
dfLine.human <- dfVal %>%
  group_by(d,r) %>%
  arrange(s) %>%
  filter(s > 1e-4) %>%
  filter(row_number() == 1);

suppl.fig6b<- dfVal %>%
  ggplot(aes(x=m,y=r,fill=log10(s))) +
  scale_y_log10( "Fold reduction of mutation rate\nby nascent RNA folding" ,position="left") +
  scale_x_log10( "Observed mutation rate \n (after reduction by nascent RNA folding)",
                 breaks = c(1e-6,1e-8,1e-10), labels=expression(10^-6,10^-8,10^-10) ) +
  geom_tile() +
  geom_line(aes(y=r,x=m),linetype=1,data=dfLine.yeast ) +
  geom_line(aes(y=r,x=m),linetype=2,data=dfLine.human ) +
  theme_classic() +
  mystyle.print() +
  scale_fill_gradientn("Selective advantage conferred\nby nascent RNA folding",
                       colours = rev(brewer.pal(11,"Spectral")),na.value=brewer.pal(11,"Spectral")[11],
                       breaks=c(-3,-7,-11),labels=expression(10^-3,10^-7,10^-11),
                       guide=guide_colorbar(title.position = "top",direction="horizontal",title.hjust=0.5,
                                            barwidth =10,ticks.colour="black",ticks.linewidth = 1 ))+
  geom_point(aes(x=x,y=y),fill="black",shape=15,data=data.frame(x=3.3e-9,y=7.26597848),size=3) +
  geom_point(aes(x=x,y=y),fill="black",shape=17,data=data.frame(x=3.3e-9,y=48.07214750),size=3) +
  theme(legend.position = "bottom")

suppl.fig6b
  
ggarrange(suppl.fig6b,suppl.fig6a,widths = c(2,1.5),labels = c("a","b"))                                          

