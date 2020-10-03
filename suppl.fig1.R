
load("mystyle.print.RData")
dfdat<-data.frame(dfmean=c(1,24),dfsd=c(0,4),dftype=c("Cytoplasm","Chromatin"),
                  dffa=factor(c(1,2)))
 suppl.fig1b<-ggplot(dfdat,aes(x=dffa,y=dfmean))+
   geom_errorbar(aes(ymin=dfmean-dfsd,ymax=dfmean+dfsd))+
   geom_bar(stat = "identity")+
 scale_x_discrete("",labels = c("Cytoplasm","Chromatin"))+
   ylab("nascent RNA/mature RNA")+
   mystyle.print()
 