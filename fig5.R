load("human_cancer/ichsape_human_ch_vivo.RData")
load("human_cancer/tcga/tcga_mut_freq.RData")

ensembl_info<-read_tsv("human_cancer/human_all_gene_ensembl.txt",col_names = T)
colnames(ensembl_info)<-c("ensg","enst","gene","gene_type","transcript_type","gene_len")

dfreact_freq<-merge(foldingStr,ensembl_info)%>%group_by(gene)%>%
  dplyr::summarise(react=-mean(react),react_gini=mean(react_gini))%>%
  merge(mut_freq)%>%filter(freq>0)

cor.test(~react+freq,dfreact_freq,method="s")
#n=428,rho=-0.13,p=0.008
cor.test(~react_gini+freq,dfreact_freq%>%filter(react_gini>0.05),method="s")
#n=428,rho=-0.14,p=0.005

###TMB
load("human_cancer/tcga/dftmb.RData")
dfreact_tmb<-merge(foldingStr,ensembl_info)%>%group_by(gene)%>%
  dplyr::summarise(react=-mean(react),react_gini=mean(react_gini))%>%
  merge(dftmb)
cor.test(~react+tmb,dfreact_tmb,method="s")
#N=689,rho=-0.18,p-value = 3.285e-06
cor.test(~react_gini+tmb,dfreact_tmb%>%filter(react_gini>0.05),method="s")
#n=689,rho=-0.20,p-value = 8.467e-08


#绘图
{
fig5a<-ggplot(data=dfreact_freq, aes(y=freq,x=react)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous(
    "Nascent RNA folding strength \n(Negative of single−stranded score)") +
  scale_y_continuous("Probability of mutated in cancer")+
  annotate("text",y=0.04,x=-0.3,label=expression("rho = -0.13  P < 0.01 "))+
  mystyle.print()

fig5b<-ggplot(data=dfreact_freq%>%filter(react_gini>0), aes(y=freq,x=react_gini)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("Nascent RNA folding strength \n(Gini index of single−stranded score)") +
  scale_y_continuous("Probability of mutated in cancer")+
  annotate("text",y=0.045,x=0.4,label=expression("rho = -0.14  P < 0.01 "))+
  mystyle.print()

fig5c<-ggplot(data=dfreact_tmb, aes(y=tmb,x=react)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("Nascent RNA folding strength \n(Negative of single−stranded score)") +
  scale_y_continuous("Mutation density in cancer")+
  annotate("text",y=0.15,x=-0.3,label=expression("rho = -0.18  P < " * 10^-{5}))+
  mystyle.print()


fig5d<-ggplot(data=dfreact_tmb%>%filter(react_gini>0.05), aes(y=tmb,x=react_gini)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("Nascent RNA folding strength \n(Gini index of single−stranded score)") +
  scale_y_continuous("Mutation density in cancer")+
  annotate("text",y=0.1,x=0.4,label=expression("rho = -0.20  P < "*10^-{7}))+
  mystyle.print()


ggarrange(fig5a,fig5b,fig5c,fig5d,ncol = 2,nrow = 2,
          labels = c("a","b","c","d"))
}
