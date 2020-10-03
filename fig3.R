library(ggpubr)

{
load("R-loop/gene.r.score.RData")
load("RData.file/modif_react_vivo_final.RData")
load("mutation/dfmutation_eachgene.RData")


allGene<-union(names(modif_react),names(gene.r.score))%>%
  union(names(dfmutation_eachgene))

react_r.score_mut<-mclapply(allGene,function(thisGene){
  mut<-dfmutation_eachgene[[thisGene]] %>%getElement("mutRate")%>%mean()
  react<- modif_react[[thisGene]]%>%getElement("RSR")%>%mean()
  react<-(-react)
  react_gini<-modif_react[[thisGene]] %>%getElement("RSR") %>%Gini()
  #r.score<-gene.r.score[[thisGene]] %>%getElement("r.score")%>%mean()
  if(is.null(gene.r.score[[thisGene]])){r.score<-NaN}else
  {r.score<-gene.r.score[[thisGene]] %>%getElement("r.score")}
  return(data.frame(gene=thisGene,mut,react, react_gini,r.score))
},mc.cores = 40)%>%rbind.fill()

cor.test(~react+mut,react_r.score_mut,method = "s",na.action = na.omit)
#N=1728, rho=-0.22 ,p<10e-19
cor.test(~react+r.score,react_r.score_mut,method = "s",na.action = na.omit)
#N=1520,  rho=-0.19,p<10e-12
cor.test(~mut+r.score,react_r.score_mut,method = "s",na.action = na.omit)
#N=5691, rho=0.21,p<10e-54

cor.test(~react_gini+mut,react_r.score_mut,method = "s",na.action = na.omit)
#N=1728, rho=-0.16 ,p<10e-10
cor.test(~react_gini+r.score,react_r.score_mut,method="s",na.action = na.omit)
#N=1520,  rho=-0.35,p<10e-45

}
####################

dat.fig4a<-react_r.score_mut%>%dplyr::select(react,r.score)%>%na.omit()
fig4a<- ggplot(dat.fig4a,aes(x=react,y=r.score)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("Nascent RNA folding strength \n (Negative of single-stranded score)") +
  scale_y_continuous(limits = c(0.8,1.12),"R-loop score")+
  annotate("text",x=-1.2,y=1.1,label=expression("rho = -0.19  P < " * 10^{-12}))+
  mystyle.print()
  
dat.fig4b<-react_r.score_mut%>%dplyr::select(react_gini,r.score)%>%na.omit()
fig4b<- ggplot(dat.fig4b,aes(x=react_gini,y=r.score)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("Nascent RNA folding strength \n (Gini index of single-stranded score)") +
  scale_y_continuous(limits = c(0.8,1.12),"R-loop score")+
  annotate("text",x=0.3,y=1.1,label=expression("rho = -0.35  P < " * 10^{-45}))+
  mystyle.print()

dat.fig4c<-react_r.score_mut%>%dplyr::select(mut,r.score)%>%na.omit()
fig4c<- ggplot(dat.fig4c,aes(x=r.score,y=mut)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("R-loop score") +
  scale_y_continuous(limits = c(0,0.95),"Mutation rate")+
  annotate("text",x=0.9,y=0.9,label=expression("rho = 0.21  P < " * 10^{-54}))+
  mystyle.print()

dat.fig4d<-react_r.score_mut%>%dplyr::select(mut,react)%>%na.omit()
fig4d<- ggplot(dat.fig4d,aes(x=react,y=mut)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("Nascent RNA folding strength \n (Negation of single-stranded score") +
  scale_y_continuous(limits = c(0,0.9),"Mutation rate")+
  annotate("text",x=-1.5,y=0.8,label=expression("rho = -0.22  P < " * 10^{-19}))+
  mystyle.print()

dat.fig4e<-react_r.score_mut%>%dplyr::select(mut,react_gini)%>%na.omit()
fig4e<- ggplot(dat.fig4e,aes(x=react_gini,y=mut)) +
  geom_point()+
  geom_smooth(method = 'lm',se=F,color='red',size=1)+
  scale_x_continuous("Nascent RNA folding strength \n (Gini index of single-stranded score)") +
  scale_y_continuous(limits = c(0,0.9),"Mutation rate")+
  annotate("text",x=0.25,y=0.8,label=expression("rho = -0.16  P < " * 10^{-10}))+
  mystyle.print()


{
  load("R-loop/base.r.score.RData")
  #每个位点进行MHtest分析
  allGene<-union(names(modif_react),names(base.r.score))%>%
    union(names(dfmutation_eachgene))
  
  ####所有基因总的OR
  {
    dfallGeneOR=data.frame()
    
    x="r.score";y="RSR"
    x="mutRate";y="r.score"
    x="mutRate";y="RSR"
    {
      df.MHtest <- mclapply(allGene, function(thisGene){
        if(x=="r.score"){df.a<-base.r.score[[thisGene]];if(is.null(df.a)){return(NULL)}}
        if(x=="mutRate"){
          df.a<-dfmutation_eachgene[[thisGene]];if(is.null(df.a)){return(NULL)}
          if(nrow(filter(df.a,mutRate>0))< 20){return(NULL)}
          df.a<-df.a%>% filter(eval(parse(text = sprintf("%s > 0",x))))
        }
        if(y=="RSR"){df.b<-modif_react[[thisGene]]
        # if(is.null(df.b)){return(NULL)}else{df.b$RSR<- -df.b$RSR}
        if(is.null(df.b)){return(NULL)}else{df.b$RSR<- df.b$RSR}
        }
        if(y=="r.score"){df.b<-base.r.score[[thisGene]]}
        if(is.null(df.b)){return(NULL)}
        
        myMerge<- merge(df.a,df.b,by="base")
        if(nrow(myMerge)<10){return(NULL)}
        median_x <- mean(eval(parse(text = sprintf("myMerge$%s",x))))
        less_x <- myMerge %>% filter(eval(parse(text = sprintf("%s <= median_x",x)))) %>% getElement("base")
        ge_x <-  myMerge %>% filter(eval(parse(text = sprintf("%s > median_x",x)))) %>% getElement("base")
        
        
        median_y <- mean(eval(parse(text = sprintf("myMerge$%s",y))))
        less_y <- myMerge %>% filter(eval(parse(text = sprintf("%s <= median_y",y))))%>%getElement("base")
        ge_y <-  myMerge %>% filter(eval(parse(text = sprintf("%s > median_y",y)))) %>%getElement("base")  
        
        
        a<- length(intersect(ge_x,ge_y))
        b <- length(intersect(ge_x,less_y))
        c <- length(intersect(less_x,ge_y))  
        d <- length(intersect(less_x,less_y))
        
        data.frame(thisGene,a,b,c,d)                                                 
      },mc.cores = 40)%>%rbind.fill
      ###
      mydata.frame <- df.MHtest[1:nrow(df.MHtest),2:5] 
      mydata.frame <- mydata.frame+1
      my.matrix <- as.matrix(mydata.frame)
      myvector <- as.vector(t(my.matrix))
      dims=c(2,2,nrow(mydata.frame));
      myarray <- array(myvector,dims)
      myOR<-mantelhaen.test(myarray)
      
      ###
      # x="mutRate";y="RSR"
      # common odds ratio= 1.082635, p= 1.245e-06
      # or>1表示碱基越趋于单链，突变率越大
      
      # x="r.score";y="RSR"
      # common odds ratio= 1.088557, p=2.336383e-84
      # or>1表示碱基越趋于单链，r.score越大
      
      # x="mutRate";y="r.score"
      # common odds ratio= 1.16287, p=2.545722e-29
      # or>1表示碱基越趋于单链，突变率越大
      
      df_1000_sd<- mclapply(1:1000, function(x){
        sample_data<- sample_n(mydata.frame,nrow(mydata.frame),replace = T);
        my.matrix <- as.matrix(sample_data)
        myvector <- as.vector(t(my.matrix))
        dims=c(2,2,nrow(mydata.frame));
        myarray <- array(myvector,dims)
        mytest<- mantelhaen.test(myarray)
        return(data.frame(commonOR=mytest$estimate))
      },mc.cores = 40)%>% rbind.fill()
      my.sd<-sd(df_1000_sd$commonOR)
      myallGeneOR=data.frame(or=myOR$estimate,sd=my.sd,type=sprintf("%s-%s",x,y))
    }
    dfallGeneOR<-rbind(myallGeneOR,dfallGeneOR)
    dfallGeneOR$name<-c("OR1",
                        "OR2",
                        "OR3")
    
    save(dfallGeneOR,fig4g,file = "RData.file/fig4g.RData")
    fig4g<-ggplot(dfallGeneOR,aes(x=name,y=or))+
      geom_errorbar(aes(ymax=or+sd,ymin=or-sd,width=0.3))+
      geom_bar(stat = "identity",position=position_dodge(),width=0.6,fill="grey80",colour = "black") +
      #scale_fill_manual(values = c("grey80","grey80"))+
      #coord_flip()+
      labs(y="Combined odds ratio",x="")+
      theme(legend.position = "")+
      #theme(axis.text.x=element_text(size = 14,colour = "black",angle=45,vjust = 0.8,hjust=0.8))+
      geom_hline(aes(yintercept=1),lty=2)+
      annotate("text",x=1,y=1.25,label=expression("P <" *" "*10^{-83}))+
      annotate("text",x=2,y=1.25,label=expression("P <" *" "*10^{-29}))+
      annotate("text",x=3,y=1.25,label=expression("P <" *" "*10^{-5}))+
      coord_cartesian(ylim = c(0.9,1.3))+
      scale_y_continuous(breaks = c(0.9,1.0,1.1,1.2))+
      mystyle.print()
    
    
    }
  
  ######################
  #分子内相关性分析
  x="mutRate";y="RSR"
  x="r.score";y="RSR"
  x="mutRate";y="r.score"
  thisGene<-allGene[167]
  {
    df.1000.set<-mclapply(1:1000, function(i){
      df.cortest <- mclapply(allGene, function(thisGene){
        if(x=="r.score"){df.a<-base.r.score[[thisGene]];if(is.null(df.a)){return(NULL)}}
        if(x=="mutRate"){df.a<-dfmutation_eachgene[[thisGene]]
        if(is.null(df.a)){return(NULL)}
        if(nrow(filter(df.a,mutRate>0))< 20){return(NULL)}
        df.a<-df.a%>% filter(eval(parse(text = sprintf("%s > 0",x))))
        }
        
        if(y=="RSR"){df.b<-modif_react[[thisGene]]}
        if(is.null(df.b)){return(NULL)}else{df.b$RSR<- -df.b$RSR}
        if(y=="r.score"){df.b<-base.r.score[[thisGene]]}
        if(is.null(df.b)){return(NULL)}
        
        myMerge<- merge(df.a,df.b,by="base")
        if(nrow(myMerge)<10){return(NULL)}
        random_df.a<-sample(eval(parse(text = sprintf("myMerge$%s",x))),
                            length(eval(parse(text = sprintf("myMerge$%s",x)))),replace = F)
        mycor<-cor.test(eval(parse(text = sprintf("myMerge$%s",x))),
                        eval(parse(text = sprintf("myMerge$%s",y))),method = "s")
        random_cor<-cor.test(eval(parse(text = sprintf("myMerge$%s",y))),
                             random_df.a ,method = "s")
        return(data.frame(cor=mycor$estimate,
                          random_cor=random_cor$estimate,gene=thisGene))
        
      },mc.cores = 25)%>%rbind.fill()
      
      return(data.frame(cor=mean(df.cortest$cor,na.rm=T),
                        random_cor=mean(df.cortest$random_cor,na.rm=T)))
    },mc.cores = 20)%>%rbind.fill()
    save(df.1000.set,file = sprintf("RData.file/df.1000.set_%s_%s.RData",x,y))
    
    
    
  }
  ####画图
  
  
  
  {
    x="r.score";y="RSR"
    load(sprintf("RData.file/df.1000.set_%s_%s.RData",x,y))
    df.1000.set$cor<--df.1000.set$cor;df.1000.set$random_cor<- -df.1000.set$random_cor
    observed_1<- mean(df.1000.set$cor)
    pvalue_1=(df.1000.set%>%filter(random_cor<observed_1)%>%nrow())/1000
    fig4h_1<- ggplot(df.1000.set,aes(x=random_cor))+geom_histogram(fill="white",colour="black")+
      scale_x_continuous("Mean correlation between negative of \nsingle-stranded score and R-loop score")+
      scale_y_continuous("Frequency")+
      coord_flip()+
      theme(plot.margin=unit(rep(1,4),'lines'))+
      annotate('segment', x=observed_1, xend=observed_1, y=100, yend=0, arrow=arrow(),size=0.5)+
      annotate('text', x=-0.008, y=50,label="P = 0 ",size=4)+
      mystyle.print()
    
    x="mutRate";y="r.score"
    load(sprintf("RData.file/df.1000.set_%s_%s.RData",x,y))
    observed_2<-mean(df.1000.set$cor)
    pvalue_2=(df.1000.set%>%filter(random_cor>observed_2)%>%nrow())/1000
    fig4h_2<- ggplot(df.1000.set,aes(x=random_cor))+geom_histogram(fill="white",colour="black")+
      scale_x_continuous("Mean correlation between R-loop\nscore and mutation rate")+
      scale_y_continuous("Frequency")+
      coord_flip()+
      theme(plot.margin=unit(rep(1,4),'lines'))+
      annotate('segment', x=observed_2, xend=observed_2, y=75, yend=0, arrow=arrow(),size=0.5)+
      annotate('text', x=0.01, y=50,label="P = 0.019 ",size=4)+
      mystyle.print()
    
    x="mutRate";y="RSR"
    load(sprintf("RData.file/df.1000.set_%s_%s.RData",x,y))
    df.1000.set$cor<--df.1000.set$cor;df.1000.set$random_cor<- -df.1000.set$random_cor
    observed_3<-mean(df.1000.set$cor)
    pvalue_3=(df.1000.set%>%filter(random_cor<observed_3)%>%nrow())/1000
    fig4h_3<- ggplot(df.1000.set,aes(x=random_cor))+geom_histogram(fill="white",colour="black")+
      scale_x_continuous("Mean correlation between negative of \nsingle-stranded score and mutation rate")+
      scale_y_continuous("Frequency")+
      coord_flip()+
      theme(plot.margin=unit(rep(1,4),'lines'))+
      annotate('segment', x=observed_3, xend=observed_3, y=70, yend=0, arrow=arrow(),size=0.5)+
      annotate('text', x=-0.01, y=25,label="P = 0.037 ",size=4)+
      mystyle.print()
  }
  
  
  
  
  
  ##########################
  x="mutRate";y="RSR"
  dfmut_RSR<-df.cortest%>%filter(cor>=quantile(df.cortest$cor,0.90)) 
  x="r.score";y="RSR"
  dfr.score_RSR<-df.cortest%>%filter(cor>=quantile(df.cortest$cor,0.90,na.rm=T))
  x="mutRate";y="r.score"
  dfr.score_mutRate<-df.cortest%>%filter(cor>=quantile(df.cortest$cor,0.90,na.rm=T))
  
  gene<-intersect(dfmut_RSR$gene,dfr.score_RSR$gene)%>%intersect(dfr.score_mutRate$gene)
  thisGene<-gene[1]
  
  RSR<-modif_react[[thisGene]]
  colnames(RSR)[2]<-"myvalue";RSR$mytype<-"icSHAPE reactivity";RSR$dffactor<-factor(1)
  RSR$myvalue<- rollapply(RSR$myvalue, width = 50,by = 1, FUN = mean, align = "center",partial = T)
  #RSR$myvalue<-SlidingWindon(RSR$myvalue,50,5)
  
  mut<-dfmutation_eachgene[[thisGene]]
  colnames(mut)[2]<-"myvalue";mut$mytype<-"Mutation rate";mut$dffactor<-factor(3)
  mut$myvalue<- rollapply(mut$myvalue, width = 50,by = 1, FUN = mean, align = "center",partial = T)
  #mut$myvalue<-SlidingWindon(mut$myvalue,50,5)
  
  r.score<-base.r.score[[thisGene]]
  dfdata.frame<-data.frame(base=1:nrow(RSR),r.score=0,gene=thisGene)
  dfdata.frame$r.score<-r.score$r.score[match(dfdata.frame$base,r.score$base)]
  dfdata.frame[is.na(dfdata.frame)] <- 0
  r.score<-dfdata.frame
  colnames(r.score)[2]<-"myvalue";r.score$mytype<-"R-loop score";r.score$dffactor<-factor(2)
  r.score$myvalue<- rollapply(r.score$myvalue, width = 50,by = 1, FUN = mean, 
                              align = "center",partial = T)
  #r.socre$myvalue<-SlidingWindon(r.socre$myvalue,50,5)
  
  dfdat<-rbind(RSR,r.score)%>%rbind(mut)
  #is.na(dfdat)<-0
  
  fig4f<-ggplot(dfdat,aes(x=base,y=myvalue))+geom_line(stat = "identity")+
    facet_grid(dffactor~.,scales = 'free')+
    mystyle.print()+
    scale_x_continuous(expand = c(0,0)," ",limits = c(1,1700),breaks = c(1,850,1700),
                       labels =c(1,850,1700))+
    theme(strip.text = element_blank(),
          legend.position = "none")
  }


############


p1<-ggarrange(fig3a,fig3b,fig3c,ncol = 3,labels = c("a","b","c"),align = "h")
p2<-ggarrange(fig3d,fig3e,fig3g,ncol = 3,labels = c("d","e","g"),align = "h")
p3<-ggarrange(fig3f,fig3h_1,fig3h_2,fig3h_3,ncol = 4,labels = c("f","h"),align = "hv")

ggarrange(p1,p2,p3,nrow = 3)




