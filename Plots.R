library(readxl)
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)


fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)




#Multinornal QQ plot

data_org=readRDS(paste(path_out,"/Data/The_Data",".rds",sep=""))
x<-data_org[,c('ANC','PLT')]

mqqnorm(x, main = "Multi-normal Q-Q Plot (Y_a,Y_b)")


data_org=readRDS(paste(path_out,"/Data/The_Data",".rds",sep=""))
d_sub=data_org%>%
  group_by(Child.Id)%>%
  summarise(Mean_ANC=mean(ANC),Mean_PLT=mean(PLT),
            Relapsed=ifelse((Relapsed[1])==0,'No','Yes'))



#Inital check
lnty=c("dotted","21","4452","41","solid")[c(1,5)]
p0<-NULL
colour_not=c("black","black")
p0<-ggplot(data = d_sub, aes(x = Mean_ANC, y = Mean_PLT)) +
  theme_bw()+theme(aspect.ratio = 1,panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())+
  scale_y_continuous(labels=c(expression(Q[b](0.1)),expression(Q[b](0.25)),expression(Q[b](0.5)),expression(Q[b](0.75)),expression(Q[b](0.9))),
                     breaks=c(quantile(d_sub$Mean_PLT, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  scale_x_continuous(labels=c(expression(Q[a](0.1)),expression(Q[a](0.25)),expression(Q[a](0.5)),expression(Q[a](0.75)),expression(Q[a](0.9))),
                     breaks=c(quantile(d_sub$Mean_ANC, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  #geom_point(aes(color=Relapsed))+ scale_color_manual(values=c("red","blue"))+
  geom_density2d(aes(linetype=Relapsed))+
  scale_linetype_manual(values=lnty)+
  labs(x='Y_a', y = 'Y_b')+#+xlim(c(6.5,8.25))+ylim(c(11,13.5))+
  guides(color= guide_legend("Density"),linetype= guide_legend("Relapsed", override.aes = list(color="black")), size=2)+
  theme(axis.title.x = element_text(size=12, face="bold", color=colour_not[1]),
        axis.title.y = element_text( size=12, face="bold",color=colour_not[2]),
        axis.text.x=element_text(size=11 , face= "bold", color=colour_not[1], angle=-22.5),
        axis.text.y=element_text(size=11 , face= "bold", color=colour_not[2]),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"))
p0                 
                 
#Contour                 
                 

lnty=c("dotted","21","4452","41","dashed","solid")[c(2,6)]
colour_not=c("brown","#22771c")

vect=c(0.05,0.1,0.15,0.2,0.225,0.3)
p4<-NULL
colour_not=c("black","black")
p4<-ggplot(data = d_sub, aes(x = Mean_ANC, y = Mean_PLT)) +
  theme_bw()+theme(aspect.ratio = 1,panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())+
  scale_y_continuous(labels=c(expression(Q[b](0.1)),expression(Q[b](0.25)),expression(Q[b](0.5)),expression(Q[b](0.75)),expression(Q[b](0.9))),
                     breaks=c(quantile(d_sub$Mean_PLT, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  scale_x_continuous(labels=c(expression(Q[a](0.1)),expression(Q[a](0.25)),expression(Q[a](0.5)),expression(Q[a](0.75)),expression(Q[a](0.9))),
                     breaks=c(quantile(d_sub$Mean_ANC, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  #geom_point(aes(color=Relapsed))+ scale_color_manual(values=c("red","blue"))+
  geom_density2d(aes(linetype=Relapsed, color=as.factor(..level..)), breaks=vect, size=1.2, alpha=1 )+ 
  scale_color_manual(values=c("red","yellow2","blue","darkgrey","maroon1","purple"))+ 
  scale_linetype_manual(values=lnty)+
  labs(x='Y_a', y = 'Y_b')+#+xlim(c(6.5,8.25))+ylim(c(11,13.5))+
  guides(color= guide_legend("Density"),linetype= guide_legend("Relapsed", override.aes = list(color="black")), size=2)+
  theme(axis.title.x = element_text(size=12, face="bold", color=colour_not[1]),
        axis.title.y = element_text( size=12, face="bold",color=colour_not[2]),
        axis.text.x=element_text(size=11 , face= "bold", color=colour_not[1], angle=-22.5),
        axis.text.y=element_text(size=11 , face= "bold", color=colour_not[2]),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        #legend.key.size = unit(0.5, 'cm'),
        legend.text=element_text(size=11),legend.key.size = unit(0.5, 'cm'),
        legend.key.height= unit(0.8, 'cm'),
        legend.key.width= unit(1, 'cm'))

q_anc=quantile(d_sub$Mean_ANC, probs=c(0.1,0.25,0.5,0.75,0.9))
q_plt=quantile(d_sub$Mean_PLT, probs=c(0.1,0.25,0.5,0.75,0.9))

p41<-NULL
p41<-p4+geom_segment(aes(x = mean(q_anc[c(2,3)]), y = mean(q_plt[c(4,5)]), xend = q_anc[1]-0.15, yend = mean(q_plt[c(4,5)])),
                     arrow = arrow(length = unit(0.3, "cm")), color="brown", size=1.1)

p42<-p41+geom_segment(aes(x = mean(q_anc[c(4,5)])-0.5, y = mean(q_plt[c(2,3)])+0.3, xend = q_anc[5]+0.7, yend = mean(q_plt[c(2,3)])+0.3),
                      arrow = arrow(length = unit(0.3, "cm")), color="green", size=1.1)

p43<-p42+geom_segment(aes(x = mean(q_anc[c(3,4)]), y = mean(q_plt[c(2,3)]), xend = mean(q_anc[c(3,4)])+0.5, yend = mean(q_plt[c(2,3)])-1),
                      arrow = arrow(length = unit(0.3, "cm")), color="orange", size=1.1)


p44<-p43+geom_segment(aes(x = mean(q_anc[c(3,4)]), y = mean(q_plt[c(2,3)])+0.1, xend = mean(q_anc[c(3,4)]), yend = mean(q_plt[c(2,3)])+1.2),
                      arrow = arrow(length = unit(0.3, "cm")), color="black", size=1.1)
p44




T_set=matrix(c(0.25,0.25,
               0.25,0.75,
               0.50,0.50,
               0.75,0.25,
               0.75,0.75),ncol=2, byrow=T)



#psi

library(readxl)
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)


fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)


prm=c("psi")
Par_ANC=NULL
Par_PLT=NULL
for ( x in 1:nrow(T_set))
{
  qr=readRDS(paste(path_out,"/Results/QRJM_",paste0(T_set[x,]*100, collapse="_"),".rds", sep=""))
  
  z_a=c(qr$q2.5[[prm]][1],qr$mean[[prm]][1],qr$q97.5[[prm]][1],paste0(T_set[x,]*100, collapse="_"))
  z_p=c(qr$q2.5[[prm]][2],qr$mean[[prm]][2],qr$q97.5[[prm]][2],paste0(T_set[x,]*100, collapse="_"))
  
  Par_ANC=rbind(Par_ANC,z_a)
  Par_PLT=rbind(Par_PLT,z_p)
  x=x+1
  
  
}

Par_ANC=data.frame(Par_ANC)
Par_PLT=data.frame(Par_PLT)
colnames(Par_ANC)=c("lci","mean","uci","tau")
colnames(Par_PLT)=c("lci","mean","uci","tau")
Par_ANC[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_ANC[,x]))
Par_PLT[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_PLT[,x]))
rownames(Par_ANC)=1:nrow(T_set)
rownames(Par_PLT)=1:nrow(T_set)





Par_ANC$Response="Y_a"
Par_PLT$Response="Y_b"
Par=data.frame(rbind(Par_ANC,Par_PLT))



p<- ggplot(Par, aes(x=tau, y=mean, group=Response, color=Response)) + scale_color_manual(values=c("red","blue"))+
  geom_line(size=1.2) +
  geom_point(size=2)+
  geom_hline(yintercept = 0, col = "black")+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=90))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold", angle = 90),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Association Estimates")+xlab("tau (Y_a_Y_b)")

p









#beta medicine
library(readxl)
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)


fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)



prm=c("b")
Par_ANC=NULL
Par_PLT=NULL


medicine=c("Med1","Med2")

for ( x in 1:nrow(T_set))
{
  qr=readRDS(paste(path_out,"/Results/QRJM_",paste0(T_set[x,]*100, collapse="_"),".rds", sep=""))
  
  z_a_mp=c(qr$q2.5[[prm]][1,3],qr$mean[[prm]][1,3],qr$q97.5[[prm]][1,3],paste0(T_set[x,]*100, collapse="_"),"Med1")
  z_p_mp=c(qr$q2.5[[prm]][2,3],qr$mean[[prm]][2,3],qr$q97.5[[prm]][2,3],paste0(T_set[x,]*100, collapse="_"),"Med1")
  
  z_a_mtx=c(qr$q2.5[[prm]][1,4],qr$mean[[prm]][1,4],qr$q97.5[[prm]][1,4],paste0(T_set[x,]*100, collapse="_"),"Med2")
  z_p_mtx=c(qr$q2.5[[prm]][2,4],qr$mean[[prm]][2,4],qr$q97.5[[prm]][2,4],paste0(T_set[x,]*100, collapse="_"),"Med2")
  
  Par_ANC=rbind(Par_ANC,z_a_mp,z_a_mtx)
  Par_PLT=rbind(Par_PLT,z_p_mp,z_p_mtx)
  x=x+1
  
  
}

Par_ANC=data.frame(Par_ANC)
Par_PLT=data.frame(Par_PLT)
colnames(Par_ANC)=c("lci","mean","uci","tau","med")
colnames(Par_PLT)=c("lci","mean","uci","tau","med")
Par_ANC[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_ANC[,x]))
Par_PLT[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_PLT[,x]))
rownames(Par_ANC)=NULL
rownames(Par_PLT)=NULL





Par_ANC$Response="Y_a"
Par_PLT$Response="Y_b"
Par=data.frame(rbind(Par_ANC,Par_PLT))



p1<- ggplot(Par[Par$med=="Med1",], aes(x=tau, y=mean, group=Response, color=Response)) + scale_color_manual(values=c("red","blue"))+
  geom_line(size=1.2) +
  geom_hline(yintercept = 0, col = "black")+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=90))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold", angle = 90),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Med1")+xlab("tau (Y_a_Y_b)")

p1

p2<- ggplot(Par[Par$med=="Med2",], aes(x=tau, y=mean, group=Response, color=Response)) + scale_color_manual(values=c("red","blue"))+
  geom_line(size=1.2) +
  geom_hline(yintercept = 0, col = "black")+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=90))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold", angle = 90),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Med2")+xlab("tau (Y_a_Y_b)")

p2

ggarrange(p1, p2,ncol=1,nrow=2,common.legend = TRUE,legend="right")










#Significance Polarity
library(readxl)
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)


fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)



si=function(x)
{
  if(min(x)<0 & max(x)>0)
  {
    z="Not"
  }
  if(max(x)<=0)
  {
    z="Negatively"
  }
  if(min(x)>=0)
  {
    z="Positively"
  }
  
  return(z)
}
si(c(1,1))

sig=function(x)
{
  qr=readRDS(paste(path_out,"/Results/QRJM_",paste0(T_set[x,]*100, collapse="_"),".rds", sep=""))
  
  
  l=qr$q2.5
  u=qr$q97.5
  
  z1=sapply(1:2, function(i) sapply(5:10, function(j) si(c(l$b[i,j],u$b[i,j])) ) )
  colnames(z1)=c("Y_a","Y_b")
  z2=sapply(1:6 , function(i) si(c(l$theta[i],u$theta[i])))
  z1=as.matrix(cbind(z1,z2))
  colnames(z1)=c("Y_a","Y_b","Surv")
  z1=data.frame(z1)
  z1$tau=rep(paste0(T_set[x,]*100, collapse="_"),6)
  z1$var=c("Var1","Var2","Var3","Var4","Var5","Var6")
  return(z1)
  
}

df=NULL
for(x in 1:5)
{
  df=rbind(df,sig(x))
  x=x+1
}
df=data.frame(df)
rownames(df)=NULL

#Remove the Status and Risk covariates
colrs=c("#6666FF","#FFFF33","#FF3333")
df=subset(df, !df$var%in%c("Risk","Status"))
p1<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+
  geom_tile(aes(fill=Y_a),color="white", size=0.75)+
  scale_fill_manual(drop=FALSE, values=colrs,name = "Significant")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("Y_a")

p2<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+
  geom_tile(aes(fill=Y_b),color="white", size=0.75)+
  scale_fill_manual(drop=FALSE, values=colrs,name = "Significant")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("Y_b")



p4<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+
  geom_tile(aes(fill=Surv),color="white", size=0.75)+
  scale_fill_manual(drop=FALSE, values=colrs,name = "Significant")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("Relapse-Time")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend<-get_legend(p4)
p1<-p1+theme(legend.position = "none")
p2<-p2+theme(legend.position = "none")
p4<-p4+theme(legend.position = "none")

grid.arrange(p1,p2,p4,legend,ncol=2,nrow=2)







#Survival
library(readxl)
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)


fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)

dataL<-readRDS(paste(path_out,"/Data/The_Data",".rds",sep=""))

Gk_points=as.matrix(cbind(c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
                            0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                            -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                            0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329),
                          c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                            0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                            0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                            0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)))

colnames(Gk_points)=c("cuts","weights")
Gk_points=data.frame(Gk_points)
Gk_points= Gk_points[with(Gk_points, order(cuts)),]
rownames(Gk_points)=1:nrow(Gk_points)
#G-K cuts
#sk
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
sk <- Gk_points$cuts
#G-K weights
#wk
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
wk <- Gk_points$weights



Surv_id=function(t,i, qr)
{
  Pr=qr$mean$Process
  Beta=qr$mean$b
  colnames(Beta)=c('Intercept','Week','Med1','Med2','Gender','NCI','MRD','CNS','Age','Height')
  rownames(Beta)=c("ANC","PLT")
  psi=qr$mean$psi
  names(psi)=c("ANC","PLT")
  gamma=qr$mean$theta
  names(gamma)=c('Gender','NCI','MRD','CNS','Age','Height')
  lambda_0=qr$mean$gam.h0
  
  
  
  t=t/max(dataL$Censor_week)
  d=subset(dataL, dataL$id==i)
  if(t<=min(d$Week)){
    survi=1
  }
  if(t>min(d$Week))
  {
    
    tpts0=sort(as.vector((t-min(d$Week))*(sk+1)/2 + min(d$Week)))
    tpts=sort(as.vector(c(min(d$Week),(d$Time_to_Event[1]-min(d$Week))*(sk+1)/2 + min(d$Week),d$Time_to_Event[1])))
    index_h=(i-1)*(length(sk)+1)+sapply(1:length(tpts0),  function(j) ifelse(tpts0[j]>=max(tpts),(length(tpts)-1),sum(sapply(2:length(tpts), function(i) ifelse(tpts0[j]>=tpts[i-1]&tpts0[j]<tpts[i], (i-1),0)))))
    Re=Pr[index_h,1:2]
    rownames(Re)=NULL
    n_s=length(tpts0)
    
    #medicine covariate
    medi=function(id, t)
    {
      d=subset(dataL, dataL$id==id)
      d=d[with(d, order( Week)),]
      z=NULL
      if(t <=max(d$Week) & t>=min(d$Week))
      {
        for(j in 1:(nrow(d)-1))
        {
          if(t>=d[j,c("Week")] & t<d[(j+1),c("Week")])
          {
            z = d[j,c("Med1","Med2")]
          }
        }
        if(t==max(d$Week))
        {
          z=d[nrow(d),c("Med1","Med2")]
        }
      }
      else{
        z=c(0,0)
      }
      return(as.numeric(z))
    }
    d_s=data.frame(Child.Id=rep(d$id[1],n_s), 
                   id=rep(i,n_s),
                   Intercept=1,
                   Week=tpts0,
                   Med1=NA,
                   Med2=NA,
                   Gender=rep(d$Gender[1],n_s),
                   NCI=rep(d$NCI[1],n_s),
                   MRD=rep(d$MRD[1],n_s),
                   CNS=rep(d$CNS[1],n_s), 
                   Age=rep(d$Age[1],n_s),
                   Height=rep(d$Height[1],n_s),
                   ANC_RE=Pr[index_h,1],
                   PLT_RE=Pr[index_h,2],
                   Censor=rep(d$Censor[1],n_s),
                   stringsAsFactors = FALSE)
    
    d_s[,c("Med1","Med2")]=t(sapply(1:nrow(d_s), function(i) medi(id=d_s$id[i],t=d_s$Week[i])))
    
    
    #Terms in hazard log hazard
    #constant
    a_i=lambda_0+as.numeric(gamma%*%t(d_s[1,names(gamma)]))
    #terms with time
    d_s$ANC=as.vector(as.matrix(d_s[,colnames(Beta)])%*%Beta['ANC',]+d_s$ANC_RE)
    d_s$PLT=as.vector(as.matrix(d_s[,colnames(Beta)])%*%Beta['PLT',]+d_s$PLT_RE)
    d_s$v= psi['ANC']*d_s$ANC +psi['PLT']*d_s$PLT
    
    z1=exp(a_i)*0.5*(t-min(d$Week)) * (as.numeric(wk%*%exp(d_s$v)))
    survi=exp(-z1) 
    
  }
  
  return(survi)
}

id_list=sort(unique(dataL$id))
time_sets=seq(25,250,25)
#time_sets=c(50,100)
T_set=matrix(c(0.25,0.25,
               0.25,0.75,
               0.50,0.50,
               0.75,0.25,
               0.75,0.75),ncol=2, byrow=T)


V=NULL
for(x in 1:nrow(T_set))
{
  
  qr=readRDS(paste(path_out,"/Results/QRJM_",paste0(T_set[x,]*100, collapse="_"),".rds", sep=""))
  z1=data.frame(t(c(0,1,1,1,1)))
  colnames(z1)=c("Week","Q1","Median","Q3","Mean")
  for(j in 1:length(time_sets))
  {
    z=NULL
    z= sapply(1:length(id_list), function(i) Surv_id(t=time_sets[j],i=id_list[i],qr))
    z_sum=NULL
    z_sum=c(time_sets[j],quantile(z,probs=c(0.25,0.50,0.75)),mean(z))
    z1=rbind(z1,z_sum)
    print(j)
    j=j+1
  }
  
  z2=data.frame(z1)
  rownames(z2)=NULL
  z2$tau=paste0(T_set[x,]*100, collapse="_")
  V=rbind(V,z2)
  print(x)
  print(Sys.time())
  x=x+1
  
}
V=data.frame(V)
saveRDS(V,paste(path_out,"/QRJM_Surv.rds", sep=""))


V=readRDS(paste(path_out,"/QRJM_Surv.rds", sep=""))
cols=brewer.pal(n = nrow(T_set), name = 'Set1')
ggplot(data = V, aes(x = Week, y = Median)) +
  geom_line(aes( color=tau),size=1.25)+theme_bw()+
  labs(fill=guide_legend(title="tau(Y_a_Y_b)"))+
  ylab("1-P(relapse)")



v1=NULL
for(x in 1:nrow(T_set))
{
  tau_st=paste0(T_set[x,]*100, collapse="_")
  d=NULL
  d=subset(V, V$tau==tau_st)
  d=subset(d, d$Week>0)
  #d=subset(d, d$Week<=250)
  d$Median=predict(loess(Median~Week, data=d))
  d$Mean=predict(loess(Mean~Week, data=d))
  d$Q1=predict(loess(Q1~Week, data=d))
  d$Q3=predict(loess(Q3~Week, data=d))
  
  v1=rbind(v1,d)
  x=x+1
}
v1=data.frame(v1)
rownames(v1)=NULL
v1=rbind(v1,subset(V,V$Week==0))
cols=brewer.pal(n = nrow(T_set), name = 'Set1')
ggplot(data = v1, aes(x = Week, y = Median)) +
  geom_line(aes( color=tau),size=1.25)+theme_bw()+
  theme(axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text( size=12, face="bold",color="black"),
        axis.text.x=element_text(size=11 , face= "bold", color="black"),
        axis.text.y=element_text(size=11 , face= "bold", color="black"),
        legend.title = element_text(colour="black", size=10, face="bold"))+labs(colour="tau(Y_a_Y_b)")+
  ylab("1-P(relapse)")







fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)

dataL<-readRDS(paste(path_out,"/Data/The_Data",".rds",sep=""))

id_list=sort(unique(dataL$id))
time_sets=seq(25,250,25)
#time_sets=c(50,100)
T_set=matrix(c(0.25,0.25,
               0.25,0.75,
               0.50,0.50,
               0.75,0.25,
               0.75,0.75),ncol=2, byrow=T)




#Longitudinal plots
Long_id=function(t,i, qr)
{
  Pr=qr$mean$Process
  Beta=qr$mean$b
  colnames(Beta)=c('Intercept','Week','Med1','Med2','Gender','NCI','MRD','CNS','Age','Height')
  rownames(Beta)=c("ANC","PLT")
  
  Gk_points=as.matrix(cbind(c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
                              0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                              -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                              0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329),
                            c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                              0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                              0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                              0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)))
  
  colnames(Gk_points)=c("cuts","weights")
  Gk_points=data.frame(Gk_points)
  Gk_points= Gk_points[with(Gk_points, order(cuts)),]
  rownames(Gk_points)=1:nrow(Gk_points)
  
  sk <- Gk_points$cuts
  wk <- Gk_points$weights
  
  t=t/max(dataL$Censor_week)
  d=subset(dataL, dataL$id==i)  
  tpts0=t
  tpts=sort(as.vector(c(min(d$Week),(d$Time_to_Event[1]-min(d$Week))*(sk+1)/2 + min(d$Week),d$Time_to_Event[1])))
  index_l=(i-1)*(length(sk)+1)+sapply(1:length(tpts0),  function(j) ifelse(tpts0[j]>=max(tpts),(length(tpts)-1),sum(sapply(2:length(tpts), function(i) ifelse(tpts0[j]>=tpts[i-1]&tpts0[j]<tpts[i], (i-1),0)))))
  Re=Pr[index_l,1:2]
  rownames(Re)=NULL
  n_s=length(tpts0)
  
  #medicine covariate
  medi=function(id, t)
  {
    d=subset(dataL, dataL$id==id)
    d=d[with(d, order( Week)),]
    z=NULL
    if(t <=max(d$Week) & t>=min(d$Week))
    {
      for(j in 1:(nrow(d)-1))
      {
        if(t>=d[j,c("Week")] & t<d[(j+1),c("Week")])
        {
          z = d[j,c("Med1","Med2")]
        }
      }
      if(t==max(d$Week))
      {
        z=d[nrow(d),c("Med1","Med2")]
      }
    }
    else{
      z=c(0,0)
    }
    return(as.numeric(z))
  }
  d_s=data.frame(Child.Id=rep(d$id[1],n_s), 
                 id=rep(i,n_s),
                 Intercept=1,
                 Week=tpts0,
                 Med1=NA,
                 Med2=NA,
                 Gender=rep(d$Gender[1],n_s),
                 NCI=rep(d$NCI[1],n_s),
                 MRD=rep(d$MRD[1],n_s),
                 CNS=rep(d$CNS[1],n_s), 
                 Age=rep(d$Age[1],n_s),
                 Height=rep(d$Height[1],n_s),
                 ANC_RE=Pr[index_l,1],
                 PLT_RE=Pr[index_l,2],
                 Censor=rep(d$Censor[1],n_s),
                 stringsAsFactors = FALSE)
  
  d_s[,c("Med1","Med2")]=t(sapply(1:nrow(d_s), function(i) medi(id=d_s$id[i],t=d_s$Week[i])))
  
  
  #Terms in hazard log hazard
  
  #terms with time
  d_s$ANC=as.vector(as.matrix(d_s[,colnames(Beta)])%*%Beta['ANC',]+d_s$ANC_RE)
  d_s$PLT=as.vector(as.matrix(d_s[,colnames(Beta)])%*%Beta['PLT',]+d_s$PLT_RE)
  
  d_s1=d_s[,c('Child.Id','id','Week','ANC','PLT')]
  return(d_s1)
}

V=NULL
for(x in 1:nrow(T_set))
{
  
  qr=readRDS(paste(path_out,"/Results/QRJM_",paste0(T_set[x,]*100, collapse="_"),".rds", sep=""))
  z1=NULL
  for(i in 1:length(id_list))
  {
    z=NULL
    z=Long_id(t=time_sets,i=id_list[i],qr)
    z1=rbind(z1,z)
    print(i)
    i=i+1
  }
  
  z2=data.frame(z1)
  colnames(z2)=c('Child.Id','id','Week','ANC','PLT')
  rownames(z2)=NULL
  z2$tau=paste0(T_set[x,]*100, collapse="_")
  V=rbind(V,z2)
  print(x)
  print(Sys.time())
  x=x+1
  
}
V=data.frame(V)
saveRDS(V,paste(path_out,"/Long_Surv.rds", sep=""))



library(dplyr)
library(readxl)
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)
library(ggplot2)

V1=V%>%
  group_by(Week, tau)%>%
  summarise(trend_ANC=mean(ANC),trend_PLT=mean(PLT))

V1$Week=V1$Week*300


cols=brewer.pal(n = nrow(T_set), name = 'Set1')
p1<-ggplot(data = V1, aes(x = Week, y = trend_ANC)) +
  geom_line(aes( color=tau),size=1.25)+theme_bw()+
  theme(axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text( size=12, face="bold",color="black"),
        axis.text.x=element_text(size=11 , face= "bold", color="black"),
        axis.text.y=element_text(size=11 , face= "bold", color="black"),
        legend.title = element_text(colour="black", size=10, face="bold"))+labs(colour="tau(Y_a_Y_b)")+
  ylab("Y_a")

p2<-ggplot(data = V1, aes(x = Week, y = trend_PLT)) +
  geom_line(aes( color=tau),size=1.25)+theme_bw()+
  theme(axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text( size=12, face="bold",color="black"),
        axis.text.x=element_text(size=11 , face= "bold", color="black"),
        axis.text.y=element_text(size=11 , face= "bold", color="black"),
        legend.title = element_text(colour="black", size=10, face="bold"))+labs(colour="tau(Y_a_Y_b)")+
  ylab("Y_b")


ggarrange(p1, p2,ncol=2,nrow=1,common.legend = TRUE,legend="bottom")



#Quantile crossing



V2=V
V2$tau_ANC=NA
V2$tau_PLT=NA
V2[,c('tau_ANC','tau_PLT')]=t(sapply(1:nrow(V2),function(x) as.numeric(strsplit(V2$tau[x],split="_")[[1]])/100))


V1_anc=V2%>%
  group_by(tau_ANC)%>%
  summarise(trend_ANC=mean(ANC))
V1_plt=V2%>%
  group_by(tau_PLT)%>%
  summarise(trend_PLT=mean(PLT))

V1$Week=V1$Week*300


p1<-ggplot(data = V1_anc, aes(x = tau_ANC, y = trend_ANC))+geom_line(linetype='dashed')+
  geom_point(size=3, col='blue')+theme_bw()+
  theme(axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text( size=12, face="bold",color="black"),
        axis.text.x=element_text(size=11 , face= "bold", color="black"),
        axis.text.y=element_text(size=11 , face= "bold", color="black"),
        legend.title = element_text(colour="black", size=10, face="bold"))+
  scale_x_continuous(breaks=c(0.25,0.50,0.75))+
  ylab("Y_a")+xlab('tau')

p2<-ggplot(data = V1_plt, aes(x = tau_PLT, y = trend_PLT))+geom_line(linetype='dashed')+
  geom_point(size=3, col='blue')+theme_bw()+
  theme(axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text( size=12, face="bold",color="black"),
        axis.text.x=element_text(size=11 , face= "bold", color="black"),
        axis.text.y=element_text(size=11 , face= "bold", color="black"),
        legend.title = element_text(colour="black", size=10, face="bold"))+
  scale_x_continuous(breaks=c(0.25,0.50,0.75))+
  ylab("Y_b")+xlab('tau')


ggarrange(p1, p2,ncol=2,nrow=1,common.legend = TRUE,legend="bottom")








#Convergence plots
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)

fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)

dataL<-readRDS(paste(path_out,"/Data/The_Data",".rds",sep=""))

id_list=sort(unique(dataL$id))
#time_sets=c(50,100)
T_set=matrix(c(0.25,0.25,
               0.25,0.75,
               0.50,0.50,
               0.75,0.25,
               0.75,0.75),ncol=2, byrow=T)



#Association Convergence
response_name=c("ANC","PLT")
nCluster=nrow(T_set)
iter=1000
association=data.frame(Iterations=rep(1:iter,nCluster*length(response_name)),
                       Response=rep(rep(response_name,each=iter),nCluster),
                       Cluster=rep(1:nCluster,each=iter*length(response_name)),
                       value=NA,
                       cumulative_value=NA)

for(i in 1:nCluster)
{
  jagsfits<-NULL
  jagsfits<-readRDS(paste(path_out,"/Results/QRJM_",paste0(T_set[i,]*100, collapse="_"),".rds", sep=""))
  association[association$Cluster==i & association$Response==response_name[1],"value"]=jagsfits$sims.list$psi[,1]
  association[association$Cluster==i & association$Response==response_name[2],"value"]=jagsfits$sims.list$psi[,2]
  association[association$Cluster==i & association$Response==response_name[1],"cumulative_value"]=cummean(jagsfits$sims.list$psi[,1])
  association[association$Cluster==i & association$Response==response_name[2],"cumulative_value"]=cummean(jagsfits$sims.list$psi[,2])
  print(i)
  i=i+1
  
}

plt_list <- vector(mode='list', length=5)
for(y in 1:5)
{
  color_theme=c("red","blue")
  means=c(sapply(c(1,2), function(x) mean(association[association$Cluster==y & association$Response==response_name[x],"value"])))
  p1<-NULL
  p2<-NULL
  p4<-NULL
  p1<-ggplot(data=association[association$Cluster==y,],aes(x=value,color=Response))+
    geom_density(size=1.1)+theme_bw()+scale_color_manual(values=color_theme)+
    xlab("")+ylab("Density")+
    geom_vline(xintercept = means
               , color=color_theme, size=1.1)+
    theme(plot.title = element_text(hjust = 0.5))
  p1<-p1+theme(legend.position="none")
  p2<-ggplot(data=association[association$Cluster==y,],aes(x=Iterations,y=value, color=Response))+
    geom_line(size=1.1)+theme_bw()+scale_color_manual(values=color_theme)+
    xlab("Iterations")+ylab("Trace")+
    geom_hline(yintercept = means,
               color=color_theme, size=1.1)+
    theme(plot.title = element_text(hjust = 0.5))
  p2<-p2+theme(legend.position="none")
  
  plt_name=paste('QRJM_',paste(T_set[y,]*100,collapse="_"),sep="")
  p4<-grid.arrange(p1,p2, ncol=2, nrow=1,widths=c(1,1),top=plt_name) 
  
  plt_list[[y]]<-p4
  
  y=y+1
  
}

p1<-ggplot(data=association[association$Cluster==1,],aes(x=value,color=Response))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=color_theme)+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(0,0)
             , color=color_theme, size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend<-get_legend(p1)

p1<-plt_list[[1]]
p2<-plt_list[[2]]
p3<-plt_list[[3]]
p4<-plt_list[[4]]
p5<-plt_list[[5]]
p6<-legend

grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2, nrow=3)

