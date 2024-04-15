library(energy)
library(simFrame)
library(ggplot2)
library(cowplot)

simulation_rounds <- c(1:500)
label <- factor(c(rep('Distance Correlation',500),
                  rep('Pearson Correlation', 500),
                  rep('Spearman Correlation',500)))

corr <- function(u1){
  #u2=1.3*u1-0.4  #corr1
  #u2=0.23*(0.55*u1+0.11)^3 #corr2
  #u2=exp(0.8*u1-2.4) #corr3
  #u2=-1.6*log(u1^2+0.38)+2.4 #corr4
  u2=(0.65*sin(0.42*u1+0.38))/(cos(-1.64*u1-0.2)) #corr5
  return(u2)
}

sim <- function(n,ncol_X,ncol_Y,N_eff_col_Y,weight_Y,noiseSD,xNArate,yNArate){
  myvalue<-NULL
  nac1 <- NAControl(NArate=xNArate)
  nac2 <- NAControl(NArate=yNArate)
  a2<-t(rep(0,ncol_Y))
  for (m in 1:50000000){
    u1<-rnorm(n,0,1)
    e1<-matrix(rnorm(ncol_X*n,0,noiseSD), nrow = n)
    a1<-t(runif(ncol_X,-2,2))
    phosite<-u1%*%a1+e1
    phosite<-as.data.frame(phosite)
    phosite<-setNA(phosite, nac1)
    
    e<-rnorm(n,0,0.01)
    u2=corr(u1)+e
    e2<-matrix(rnorm(ncol_Y*n,0,noiseSD),nrow = n)
    eff_col_Y<-as.numeric(sample(1:ncol_Y,N_eff_col_Y,replace = F))
    a2[,eff_col_Y]<-weight_Y
    protein_data<-u2%*%a2+e2
    protein_data<-as.data.frame(protein_data)
    protein_data<-setNA(protein_data, nac2)
    
    array_phosite<-rowMeans(phosite,na.rm = T)
    array_protein_data<-rowMeans(protein_data,na.rm = T)
    
    Dist_phosite<-dist(phosite)
    Dist_protein_data<-dist(protein_data)
    
    Dist_phosite<-as.matrix(Dist_phosite)
    Dist_protein_data<-as.matrix(Dist_protein_data)
    
    if (!(all(is.finite(c(Dist_phosite,Dist_protein_data))))){
      Dp<-NA
      DCor<-NA
      Pp<-NA
      PCor<-NA
      Sp<-NA
      SCor<-NA
    }
    
    else {
      Dp<-dcorT.test(phosite, protein_data)$p.value
      DCor<-dcor(phosite, protein_data)
      
      Pp<-cor.test(array_phosite,array_protein_data,method="pearson")$p.value
      PCor<-cor(array_phosite,array_protein_data,method="pearson")
      
      Sp<-cor.test(array_phosite,array_protein_data,method="spearman")$p.value
      SCor<-cor(array_phosite,array_protein_data,method="spearman")
    }
    
    value<-data.frame(DCor,Dp,PCor,Pp,SCor,Sp)
    myvalue<-rbind(value,myvalue)
    rownames(myvalue)<-NULL
    myvalue<-myvalue[complete.cases(myvalue),]
    myvalue<<-myvalue
    
    if (m%%1000==0 | nrow(myvalue)==500){
      print(paste(m," rounds finished, ",nrow(myvalue),"/500 effective results obtained",sep = ""))
    }
    
    if (nrow(myvalue) >= 500)
      break
  }
  
  myvalue$DBH<-p.adjust(myvalue$Dp,method = 'BH')
  myvalue$PBH<-p.adjust(myvalue$Pp,method = 'BH')
  myvalue$SBH<-p.adjust(myvalue$Sp,method = 'BH')
  return(myvalue)
} 

result <- sim(n=50,ncol_X=50,ncol_Y=70,N_eff_col_Y=0,
              weight_Y=0,noiseSD=0.1,xNArate=0,yNArate=0)

###############################################
DATA<-result
{
  p1 <- ggplot(DATA,aes(x=DCor,color=DBH))+
    theme_classic()+
    scale_color_gradient(limits=c(0,1),low="midnightblue",high="lightskyblue")+
    geom_rug(length=unit(0.1,"npc"),sides="b")+
    guides(color=guide_legend(override.aes=list(size=7)))+
    #geom_density(color="dodgerblue3",size=0.8,fill="dodgerblue2",alpha=0.4)+
    geom_density(color="white",size=0.01,fill="gray60",alpha=0.4)+
    xlim(-1,1)+
    xlab("distance coefficients")+
    ylim(-2,30)+
    ylab("density")+
    theme(legend.title=element_text(size=16,hjust=0,vjust=0.5))+
    theme(legend.text=element_text(size=14,hjust=0,vjust=0.5))+
    theme(axis.title.x=element_text(size=20,color="black",face="plain",vjust=0.5,hjust=0.5))+
    theme(axis.text.x=element_text(size=18,color="gray50",face="bold",vjust=0.5,hjust=0.5))+
    theme(axis.title.y=element_text(size=20,color="black",face="plain",vjust=0.5,hjust=0.5))+
    theme(axis.text.y=element_text(size=18,color="gray50",face="bold",vjust=0.5,hjust=0.5))
  
  p2 <- ggplot(DATA,aes(x=PCor,color=PBH))+
    theme_classic()+
    scale_color_gradient(limits=c(0,1),low="midnightblue",high="lightskyblue")+
    geom_rug(length=unit(0.1,"npc"),sides="b")+
    guides(color=guide_legend(override.aes=list(size=7)))+
    #geom_density(color="royalblue3",size=0.8,fill="royalblue2",alpha=0.4)+
    geom_density(color="white",size=0.01,fill="gray60",alpha=0.4)+
    xlim(-1,1)+
    xlab("Pearson coefficients")+
    ylim(-0.28,4)+
    ylab("density")+
    theme(legend.title=element_text(size=16,hjust=0,vjust=0.5))+
    theme(legend.text=element_text(size=14,hjust=0,vjust=0.5))+
    theme(axis.title.x=element_text(size=20,color="black",face="plain",vjust=0.5,hjust=0.5))+
    theme(axis.text.x=element_text(size=18,color="gray50",face="bold",vjust=0.5,hjust=0.5))+
    theme(axis.title.y=element_text(size=20,color="black",face="plain",vjust=0.5,hjust=0.5))+
    theme(axis.text.y=element_text(size=18,color="gray50",face="bold",vjust=0.5,hjust=0.5))
  
  p3 <- ggplot(DATA,aes(x=SCor,color=SBH))+
    theme_classic()+
    scale_color_gradient(limits=c(0,1),low="midnightblue",high="lightskyblue")+
    geom_rug(length=unit(0.1,"npc"),sides = "b")+
    guides(color=guide_legend(override.aes=list(size=7)))+
    #geom_density(color="navy",size=0.8,fill="blue3",alpha=0.3)+
    geom_density(color="white",size=0.01,fill="gray60",alpha=0.4)+
    xlim(-1,1)+
    xlab("Spearman coefficients")+
    ylim(-0.28,4)+
    ylab("density")+
    theme(legend.title=element_text(size=16,hjust=0,vjust=0.5))+
    theme(legend.text=element_text(size=14,hjust=0,vjust=0.5))+
    theme(axis.title.x=element_text(size=20,color="black",face="plain",vjust=0.5,hjust=0.5))+
    theme(axis.text.x=element_text(size=18,color="gray50",face="bold",vjust=0.5,hjust=0.5))+
    theme(axis.title.y=element_text(size=20,color="black",face="plain",vjust=0.5,hjust=0.5))+
    theme(axis.text.y=element_text(size=18,color="gray50",face="bold",vjust=0.5,hjust=0.5))
  
  p4 <- cowplot::plot_grid(p1, p2, p3, nrow=3,labels = LETTERS[1:3])
  p4
}

{
  DCor <- abs(c(result_14_9_corr1$DCor,
                result_14_18_corr1$DCor,
                result_14_27_corr1$DCor,
                result_14_36_corr1$DCor,
                result_14_45_corr1$DCor,
                result_14_54_corr1$DCor,
                result_14_9_corr5$DCor,
                result_14_18_corr5$DCor,
                result_14_27_corr5$DCor,
                result_14_36_corr5$DCor,
                result_14_45_corr5$DCor,
                result_14_54_corr5$DCor))
  
  PCor <- abs(c(result_14_9_corr1$PCor,
                result_14_18_corr1$PCor,
                result_14_27_corr1$PCor,
                result_14_36_corr1$PCor,
                result_14_45_corr1$PCor,
                result_14_54_corr1$PCor,
                result_14_9_corr5$PCor,
                result_14_18_corr5$PCor,
                result_14_27_corr5$PCor,
                result_14_36_corr5$PCor,
                result_14_45_corr5$PCor,
                result_14_54_corr5$PCor))
  
  SCor <- abs(c(result_14_9_corr1$SCor,
                result_14_18_corr1$SCor,
                result_14_27_corr1$SCor,
                result_14_36_corr1$SCor,
                result_14_45_corr1$SCor,
                result_14_54_corr1$SCor,
                result_14_9_corr5$SCor,
                result_14_18_corr5$SCor,
                result_14_27_corr5$SCor,
                result_14_36_corr5$SCor,
                result_14_45_corr5$SCor,
                result_14_54_corr5$SCor))
  
  DBH <- c(result_14_9_corr1$DBH,
           result_14_18_corr1$DBH,
           result_14_27_corr1$DBH,
           result_14_36_corr1$DBH,
           result_14_45_corr1$DBH,
           result_14_54_corr1$DBH,
           result_14_9_corr5$DBH,
           result_14_18_corr5$DBH,
           result_14_27_corr5$DBH,
           result_14_36_corr5$DBH,
           result_14_45_corr5$DBH,
           result_14_54_corr5$DBH)

PBH <- c(result_14_9_corr1$PBH,
         result_14_18_corr1$PBH,
         result_14_27_corr1$PBH,
         result_14_36_corr1$PBH,
         result_14_45_corr1$PBH,
         result_14_54_corr1$PBH,
         result_14_9_corr5$PBH,
         result_14_18_corr5$PBH,
         result_14_27_corr5$PBH,
         result_14_36_corr5$PBH,
         result_14_45_corr5$PBH,
         result_14_54_corr5$PBH)

SBH <- c(result_14_9_corr1$SBH,
         result_14_18_corr1$SBH,
         result_14_27_corr1$SBH,
         result_14_36_corr1$SBH,
         result_14_45_corr1$SBH,
         result_14_54_corr1$SBH,
         result_14_9_corr5$SBH,
         result_14_18_corr5$SBH,
         result_14_27_corr5$SBH,
         result_14_36_corr5$SBH,
         result_14_45_corr5$SBH,
         result_14_54_corr5$SBH)

Cor <- c(DCor, PCor, SCor)
BH <- c(DBH, PBH, SBH)
}

Model <- rep(c(rep("y=9", 500), rep("y=18", 500), rep("y=27", 500),
               rep("y=36", 500), rep("y=45", 500), rep("y=54", 500)), 6)
Model <- factor(Model,
                levels = c("y=9","y=18","y=27",
                           "y=36","y=45","y=54"))

N <- rep(c(rep("Corr1", 3000), rep("Corr5", 3000)),3)
N <- factor(N,
            levels = c("Corr1",
                       "Corr5"))

Method <- c(rep("Distance", 6000),
                rep("Pearson", 6000),
                rep("Spearman", 6000))

data <- data.frame(Cor, BH, N, Model, Method)

main <- ggplot(aes(x = Cor, y = Model), data = data) +
  geom_density_ridges(aes(point_color = BH, fill = Method, alpha = Method),
                      jittered_points = TRUE, point_alpha = .5, 
                      point_size = 0.6, point_shape = 16,
                      scale = .9, size = .001, color = NA,
                      position = "raincloud", show.legend = F) +
  theme_bw() +
  scale_fill_manual(values = c("#EF0000","#384793","#008C43"))+
  scale_alpha_manual(values = c(.9,.8,.7))+
  scale_point_color_gradient(low = "gray10", high = "gray90") +
  scale_x_continuous(limits = c(0,1.05), breaks = seq(.1,.9,.2), 
                     expand = c(0,0), position = "top") +
  theme(axis.title = element_blank(),
        strip.placement = "outside",
        strip.background.x = element_blank(),
        strip.background.y = element_rect(color = "#000000", fill = "#FFFFFF"),
        strip.text.x = element_textbox(size = 12, color = "white", face = "bold",
                                       fill = "#008280", box.color = NA,
                                       halign = 0.5, linetype = 1,
                                       r = unit(5, "pt"), width = unit(1, "npc"),
                                       padding = margin(2, 0, 1, 0),
                                       margin = margin(3, 3, 3, 3)),
        strip.text.y = element_text(face = "bold", colour = "#008280", size = 16),
        axis.text.x = element_text(vjust = 5, size = 10),
        axis.text.y = element_text(size = 11, face = "bold",color = "#3F3F3F"),
        panel.grid.minor = element_blank())+
  facet_grid(Method ~ N)

{
U1<-runif(10000,-10,10)
U2<-corr(U1)
figure <- data.frame(U1,U2)
ggplot(figure, aes(x=U1,y=U2))+
  ggtitle("????量?????怨?系图示")+
  theme_bw()+
  geom_point(size=2,shape=16,color="midnightblue")+
  theme(plot.title=element_text(size=22,face="bold",hjust = 0.5))+
  theme(plot.subtitle=element_text(size=18,color='gray50',hjust = 0.5))+
  scale_y_continuous(expand=c(0,0),limits = c(-10,10))+ #????y??????
  scale_x_continuous(expand=c(0,0),limits = c(-14,28))+
  theme(axis.text.y=element_text(size=15,hjust=1,vjust=0.5,face="bold"))+
  theme(axis.text.x=element_text(size=15,hjust=0.5,vjust=0.5,face="bold"))+
  theme(axis.title.y=element_text(size=17,hjust=0.5,vjust=0.5,face="plain"))+
  theme(axis.title.x=element_text(size=17,hjust=0.5,vjust=0.5,face="plain"))+
  xlab("????量u1")+
  ylab("????量u2")
}
{
ggplot(data=result_50_nNA_lN_21_54_corr5)+
  theme_bw()+
  geom_point(aes(simulation_rounds,DCor,color=DBH,shape="Distance Correlation"),size=3)+
  geom_point(aes(simulation_rounds,PCor,color=PBH,shape="Pearson Correlation"),size=2)+
  geom_point(aes(simulation_rounds,SCor,color=SBH,shape="Spearman Correlation"),size=2)+
  scale_color_gradient(limits=c(0,1),low="midnightblue",high="lightskyblue")+
  guides(shape=guide_legend(override.aes=list(size=7)))+
  guides(color=guide_legend(override.aes=list(size=7)))+
  ggtitle("??同??????模?偷?模??实???????员?散??图",subtitle="n=50        corr5")+
  theme(plot.title=element_text(size=30,face="bold",hjust = 0.5))+
  theme(plot.subtitle=element_text(size=26,color='gray50',hjust = 0.5))+
  theme(legend.title=element_text(size=20,hjust=1,vjust=0.8))+
  theme(legend.text=element_text(size=18,hjust=0.5,vjust=0.5))+
  theme(legend.position="bottom")+
  labs(shape="??????模??",color="BH p值",label.position='top')+ 
  xlab("模??实??????")+
  theme(axis.title.x=element_text(size=28,color="black",face="bold",vjust=0.5,hjust=0.5))+
  theme(axis.text.x=element_text(size=24,color="gray50",face="bold",vjust=0.5,hjust=0.5))+
  ylab("????系??")+
  ylim(-1,1)+
  theme(axis.title.y=element_text(size=28,color="black",face="bold",vjust=0.5,hjust=0.5))+
  theme(axis.text.y=element_text(size=24,color="gray50",face="bold",vjust=0.5,hjust=0.5))
}
###############################################
Cor <- c(result_50_NA_control$DCor,result_50_NA_control$PCor,result_50_NA_control$SCor)
Cor_A <- abs(Cor)
BH <- c(result_50_NA_control$DBH,result_50_NA_control$PBH,result_50_NA_control$SBH)
data <- data.frame(simulation_rounds,label,Cor,BH)

ggplot(data,aes(label,Cor,color=BH))+
  geom_point(size=3)+
  scale_color_gradient(low = "midnightblue",high = "lightskyblue")+
  ggtitle("??同??????模?偷?模??实???????员?散??图-2")+
  theme(plot.title=element_text(size=15,face="bold",hjust = 0.5))+
  theme(plot.subtitle=element_text(size=13,color='gray50',hjust = 0.5))+
  theme(legend.position="right")+
  labs(color="BH p值",label.position='top')+ 
  xlab("??????模??")+
  theme(axis.title.x=element_text(size=17,color="black",face="plain",vjust=0.5,hjust=0.5))+
  theme(axis.text.x=element_text(size=15,color="gray50",face="bold",vjust=0.5,hjust=0.5))+
  ylab("????系??")+
  ylim(-1,1)+
  theme(axis.title.y=element_text(size=17,color="black",face="plain",vjust=0.5,hjust=0.5))+
  theme(axis.text.y=element_text(size=15,color="gray50",face="bold",vjust=0.5,hjust=0.5))

###############################################
ggplot(data,aes(label,Cor)) + 
  geom_violin(aes(fill=label),color='white',show.legend = F) +
  scale_fill_manual(values=c("lightskyblue2","steelblue3","royalblue4"))+
  geom_boxplot(width = 0.03)+
  ggtitle("??同模??????系???植?????")+
  theme(plot.title=element_text(size=15,face="bold",hjust = 0.5))+
  theme(plot.subtitle=element_text(size=13,color='gray50',hjust = 0.5))+
  xlab("??????模??")+
  theme(axis.title.x=element_text(size=17,color="black",face="plain",vjust=0.5,hjust=0.5))+
  theme(axis.text.x=element_text(size=12,color="gray50",face="bold",vjust=0.5,hjust=0.5))+
  ylab("????系??")+
  ylim(-1,1)+
  theme(axis.title.y=element_text(size=17,color="black",face="plain",vjust=0.5,hjust=0.5))+
  theme(axis.text.y=element_text(size=15,color="gray50",face="bold",vjust=0.5,hjust=0.5))

############ult$DCor,ylim=c(-1,1),pch=16)
  points(result$PCor,col='red',pch=16)
  points(result$SCor,col='brown',pch=17)
  
  plot(result$DCor,ylim=c(-1,1),pch=16)
  points(abs(result$PCor),col='red',pch=16)
  points(abs(result$SCor),col='brown',pch=17)
  
  plot(result$DBH,ylim=c(0,1),pch=16)
  points(result$PBH,col='red',pch=16)
  points(result$SBH,col='brown',pch=17)
}