#Library
library(rjags)
library(R2jags)
library(jagsUI)
library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)
library(sampling)
library(truncnorm)


fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path_out=paste(fp,"/Sim_QRJM",sep="")
setwd(path_out)
#G-K points
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

#True parameters
#---------------------------------------------------------------------------------------#
TP_thetas=function(N,dist.string_ANC,dist.string_PLT)
{
  #True Parameters
  sig1=sqrt(c(0.4,0.75))
  
  
  err_sigma=c(0.4,0.4)
  names(err_sigma)=c('ANC','PLT')
  
  Beta=NULL
  psi=NULL
  gamma=NULL
  if(dist.string_ANC=='ALD0.25' & dist.string_PLT=='ALD0.25') #Better by low ANC
  {
    Beta=matrix(c(c(0.11,0.418,      -0.295,0.027,     0.135,0.042,0.06,-0.087,0.104,0.04),
                  c(-0.907,-0.066,    0.052,0.077,    0.061,0.176,0.078,-0.008,0.078,0.041)),byrow=T, nrow=2)
    Beta=Beta*2
    psi=c(1.53,-2.62)
    gamma=c(0.752,0.384,-0.282,0.116,0.23,0.102)*2
    
    rho=c(-0.5)
    Sigma=diag(sig1[1:2])%*%matrix(c(1,rho[1],rho[1],1),byrow=T,nrow=2)%*%diag(sig1[1:2])
    colnames(Sigma)=rownames(Sigma)=c("ANC","PLT")
    
  }
  if(dist.string_ANC=='ALD0.25' & dist.string_PLT=='ALD0.75') #Best
  {
    Beta=matrix(c(c(0.073,0.35,  -0.328,0.023,     0.137,0.102,-0.019,-0.015,0.133,0.038),
                  c(0.451,-0.275,     -0.02,0.045,   -0.02,0.072,-0.041,0.055,0.041,0.01)),byrow=T, nrow=2)
    Beta=Beta*2
    psi=c(1.9,-2.3)
    gamma=c(0.743,0.968,-0.504,0.075,0.238,0.069)*2
    
    rho=c(0.8)
    Sigma=diag(sig1[1:2])%*%matrix(c(1,rho[1],rho[1],1),byrow=T,nrow=2)%*%diag(sig1[1:2])
    colnames(Sigma)=rownames(Sigma)=c("ANC","PLT")
    
  }
  
  if(dist.string_ANC=='ALD0.5' & dist.string_PLT=='ALD0.5') #Median
  {
    Beta=matrix(c(c(0.704,-0.368,  -0.285,-0.022,    0.086,0.067,0.001,-0.1,0.115,0.047),
                  c(-0.196,-0.179,  -0.005,0.086,   -0.018,0.077,0.001,0.056,0.036,0.026)),byrow=T, nrow=2)
    Beta=Beta*2
    
    psi=c(2.32,-3.13)
    gamma=c(0.625,0.536,-0.312,0.235,0.189,0.041)*2
    
    rho=c(0.5)
    Sigma=diag(sig1[1:2])%*%matrix(c(1,rho[1],rho[1],1),byrow=T,nrow=2)%*%diag(sig1[1:2])
    colnames(Sigma)=rownames(Sigma)=c("ANC","PLT")
    
  }
  if(dist.string_ANC=='ALD0.75' & dist.string_PLT=='ALD0.25') #Worst
  {
    Beta=matrix(c(c(0.0561,-0.954,   -0.243,-0.021,  0.087,0.164,-0.023,-0.035,0.126,0.063),
                  c(-0.536,0.163,    -0.037,0.076,   0.1,-0.022,0.035,0.065,-0.019,-0.027)),byrow=T, nrow=2)
    Beta=Beta*2
    psi=c(1.52,-1.96)
    gamma=c(0.791,0.172,-0.228,0.097,0.137,0.018)*2
    
    rho=c(-0.5)
    Sigma=diag(sig1[1:2])%*%matrix(c(1,rho[1],rho[1],1),byrow=T,nrow=2)%*%diag(sig1[1:2])
    colnames(Sigma)=rownames(Sigma)=c("ANC","PLT")
  }
  if(dist.string_ANC=='ALD0.75' & dist.string_PLT=='ALD0.75') #Better by high plt
  {
    Beta=matrix(c(c(0.999,-0.958,  -0.23,-0.034,   0.105,0.117,-0.002,-0.173,0.117,0.047),
                  c(0.723,-0.411,  -0.027,0.046,   -0.017,0.048,-0.033,-0.028,0.035,0.004)),byrow=T, nrow=2)
    Beta=Beta*2
    psi=c(2.24,-2.43)
    gamma=c(0.702,0.697,-0.445,0.238,0.175,0.007)*2
    
    rho=c(0.4)
    Sigma=diag(sig1[1:2])%*%matrix(c(1,rho[1],rho[1],1),byrow=T,nrow=2)%*%diag(sig1[1:2])
    colnames(Sigma)=rownames(Sigma)=c("ANC","PLT")
  }
  
  colnames(Beta)=c('Intercept','Week','Med1','Med2','Gender','NCI','MRD','CNS','Age','Height')
  rownames(Beta)=c('ANC','PLT')
  names(psi)=c('ANC','PLT')
  names(gamma)=c('Gender','NCI','MRD','CNS','Age','Height')
  lambda_0=0.05
  #Error distribution names choices are gamma, t, Normal ALD0.25,ALD0.5,ALD0.75,Unif
  #dist.string="Unif"
  
  
  #Number of id
  #N<-100
  
  #BM parts
  RE=mvrnorm(n=N*20,mu=c(0,0), Sigma=Sigma)
  RE=data.frame(RE)
  RE$id=rep(1:N, each=20)
  
  #True parameters
  Thetas<-list('N'=N,'Beta'=Beta,'Sigma'=Sigma,'RE'=RE,'err_sigma'=err_sigma,'dist.string_ANC'=dist.string_ANC,'dist.string_PLT'=dist.string_PLT,
               'psi'=psi,'gamma'=gamma,'lambda_0'=lambda_0)
  return(Thetas)
}

#Data Generation function
Data_Gen=function(Thetas)
{
  #Parameter unpacking
  Beta=Thetas$Beta
  Sigma=Thetas$Sigma
  RE=Thetas$RE
  err_sigma=Thetas$err_sigma
  psi=Thetas$psi
  gamma=Thetas$gamma
  lambda_0=Thetas$lambda_0
  dist.string_ANC=Thetas$dist.string_ANC
  dist.string_PLT=Thetas$dist.string_PLT
  N=Thetas$N
  pro_granular=nrow(RE)/N
  #id names
  id=paste("UPN_",1:N, sep="")
  
  #Number of Measurements
  N_i=sapply(1:N, function(x) as.numeric(ifelse(runif(1,0,1)<=0.3,rpois(1,3),rpois(1,10)))) +5

  #Starting times of Therapy
  t_i1=rep(1,N)
  
  #Censoring week for each id
  Censor= c((rbinom((N-1),5,0.5)*4)+280 , 300)

  #Some fixed and time-invariant covariates
  
  #Discrete
  Z1=rbinom(N,1,(0.7))
  Z2=rbinom(N,2,(0.5))*0.5
  Z3=rbinom(N,3,(0.3))*0.5
  Z4=rbinom(N,4,(0.2))*0.5
  
  #Continuous
  Z5=rtruncnorm(N, a=50, b=100, mean =75, sd = 1)
  Z5=(Z5-mean(Z5))/sd(Z5)
  Z6=rtruncnorm(N, a=1, b=19, mean =6.17, sd = 4)
  Z6=(Z6-mean(Z6))/sd(Z6)
  
  indi=function(x,a,b)
  {
    z=ifelse(x>=a & x<b,1,0)
    return(z)
  }
  
  #Process offset
  offset_Pr=(0:(N))*((nrow(RE)/N))+1
  dataL=NULL
  ProL=NULL
  ProS=NULL
  ProH=NULL
  #Process
  Pr=data.frame(matrix(rep(NA,(nrow(RE)*2)), nrow=nrow(RE)))
  for(i in 1:length(id))
  {
    #Measurement times for id[i] 
    Week=NULL
    Week=t_i1[i]+c(0,cumsum(srswr(2,(N_i[i]-1))+4))  

    d=NULL
    d=data.frame(Child.Id=rep(id[i],N_i[i]), 
                 id=rep(i,N_i[i]),
                 Intercept=1,
                 Week=Week,
                 Med1=rbeta(N_i[i], shape1=2, shape2=2, ncp = 0),
                 Med2=rbeta(N_i[i], shape1=4, shape2=2, ncp = 0),
                 Gender=rep(Z1[i],N_i[i]), 
                 NCI=rep(Z2[i],N_i[i]),
                 MRD=rep(Z3[i],N_i[i]),
                 CNS=rep(Z4[i],N_i[i]), 
                 Age=rep(Z5[i],N_i[i]),
                 Height=rep(Z6[i],N_i[i]),
                 
                 ANC_RE=NA,
                 PLT_RE=NA,
                 Censor=rep(Censor[i],N_i[i]),
                 stringsAsFactors = FALSE)
    
    tpts=sort(as.vector((d$Censor[1]-min(d$Week))*(seq(-1,1,2/pro_granular)+1)/2)) + min(d$Week)
    tpts0=sort(as.vector((d$Censor[1]-min(d$Week))*(sk+1)/2 + min(d$Week)))
    d=d[with(d, order(Week)),]
    
    d$index=(i-1)*(pro_granular)+sapply(1:nrow(d),  function(j) ifelse(d$Week[j]>=max(tpts),(length(tpts)-1),sum(sapply(2:length(tpts), function(i) ifelse(d$Week[j]>=tpts[i-1]&d$Week[j]<tpts[i], (i-1),0)))))
    index_h=(i-1)*(pro_granular)+sapply(1:length(tpts0),  function(j) ifelse(tpts0[j]>=max(tpts),(length(tpts)-1),sum(sapply(2:length(tpts), function(i) ifelse(tpts0[j]>=tpts[i-1]&tpts0[j]<tpts[i], (i-1),0)))))
    
    ProL=rbind(ProL,d[,c("id","index")])
    ProH=rbind(ProH,data.frame("id"=rep(id[i], length(index_h)),"index"=index_h))
    tpts1=sqrt(diff(tpts/max(Censor)))
    ProS =rbind(ProS,cbind(rep(d$id[1],length(tpts1)),c(1:length(tpts1))+(i-1)*(length(tpts1)),tpts1))
    colnames(ProS)=c("id","Process_time_id","Process_time")
    
    d=d[,!colnames(d)%in%c("index")]
    dataL=rbind(dataL,d)
    
    #BM Process
    #----------
    
    Pr[offset_Pr[i],] = RE[offset_Pr[i],1:2] * ProS[offset_Pr[i],3]
    
    for( w in (offset_Pr[i]+1):(offset_Pr[i+1]-1))
    {
      
      Pr[w,] = Pr[(w-1),] + RE[w,1:2] * ProS[w,3]
      w=w+1
    }
    
    i=i+1
  }
  dataL=data.frame(dataL)
  ProL=data.frame(ProL)
  ProS= data.frame(ProS)
  ProH= data.frame(ProH)
  Pr=data.frame(Pr)
  rownames(dataL)=1:nrow(dataL)
  rownames(ProL)=1:nrow(ProL)
  rownames(ProS)=1:nrow(ProS)
  rownames(ProH)=1:nrow(ProH)
  
  dataL$ANC_RE=Pr[ProL[,c("index")],1]
  dataL$PLT_RE=Pr[ProL[,c("index")],2]
  
  dataL$Intercept=1
  
  #Rearrange columns
  dataL=dataL[,c(1:2,ncol(dataL),3:(ncol(dataL)-1))]
  dataL$Week=dataL$Week/max(Censor)
  dataL$Censor=dataL$Censor/max(Censor)
  

    
    
  
  
  #Response
  
  #Error function
  Err_func=function(Num, type, scale , loc)
  {
    ald_p=function(p)
    {
      th=(1-2*p)/(p*(1-p))
      ta_sq=2/(p*(1-p))
      return(c(th, ta_sq))
    }
    
    v=NULL
    v=ald_p(as.numeric(strsplit(type,"ALD")[[1]][2]))
    #Exponential
    ex=rexp(Num, rate=1)
    z=(v[1]*ex + sqrt(v[2]*ex)*rnorm(Num,mean=0,sd=1))*scale +loc
    
    return(z)
  }
  
  #Fixed part                               #RE                                     #Error
  dataL$ANC=as.vector(as.matrix(dataL[,colnames(Beta)])%*%Beta['ANC',]+dataL$ANC_RE + Err_func(nrow(dataL),dist.string_ANC,err_sigma['ANC'],0))
  dataL$PLT=as.vector(as.matrix(dataL[,colnames(Beta)])%*%Beta['PLT',]+dataL$PLT_RE + Err_func(nrow(dataL),dist.string_PLT,err_sigma['PLT'],0))
  
  
  #Survival
  
  
  Surv_id=function(t,i)
  {
    t=t/max(Censor)
    d=subset(dataL, dataL$id==i)
    if(t<=min(d$Week)){
      survi=1
    }
    if(t>min(d$Week))
    {
      
      tpts0=sort(as.vector((t-min(d$Week))*(sk+1)/2 + min(d$Week)))
      tpts=sort(as.vector((d$Censor[1]-min(d$Week))*(seq(-1,1,2/pro_granular)+1)/2)) + min(d$Week)
      index_h=(i-1)*(pro_granular)+sapply(1:length(tpts0),  function(j) ifelse(tpts0[j]>=max(tpts),(length(tpts)-1),sum(sapply(2:length(tpts), function(i) ifelse(tpts0[j]>=tpts[i-1]&tpts0[j]<tpts[i], (i-1),0)))))
      Re=Pr[index_h,1:2]
      rownames(Re)=NULL
      n_s=length(tpts0)
      
      tpts_m=c(d$Week,(max(d$Week)+6/max(Censor)))
      index_m=sapply(1:length(tpts0),  function(j) ifelse(tpts0[j]>max(tpts_m),0,sum(sapply(2:length(tpts_m), function(i) ifelse(tpts0[j]>=tpts_m[i-1]&tpts0[j]<tpts_m[i], (i-1),0)))))
      
      d_s=data.frame(Child.Id=rep(id[i],n_s), 
                   id=rep(i,n_s),
                   Intercept=1,
                   Week=tpts0,
                   Med1=sapply(1:length(index_m),function(i1) ifelse(index_m[i1]==0, 0 ,d$Med1[index_m[i1]])),
                   Med2=sapply(1:length(index_m),function(i1) ifelse(index_m[i1]==0, 0 ,d$Med2[index_m[i1]])),
                   Gender=rep(Z1[i],n_s),
                   NCI=rep(Z2[i],n_s),
                   MRD=rep(Z3[i],n_s),
                   CNS=rep(Z4[i],n_s), 
                   Age=rep(Z5[i],n_s),
                   Height=rep(Z6[i],n_s),
                   ANC_RE=Pr[index_h,1],
                   PLT_RE=Pr[index_h,2],
                   Censor=rep(Censor[i],n_s),
                   stringsAsFactors = FALSE)
      
      
      #Terms in hazard log hazard
      #constant
      a_i=log(lambda_0)+as.numeric(gamma%*%t(d_s[1,names(gamma)]))
      #terms with time
      d_s$ANC=as.vector(as.matrix(d_s[,colnames(Beta)])%*%Beta['ANC',]+d_s$ANC_RE)
      d_s$PLT=as.vector(as.matrix(d_s[,colnames(Beta)])%*%Beta['PLT',]+d_s$PLT_RE)
      d_s$v= psi['ANC']*d_s$ANC +psi['PLT']*d_s$PLT
      
      z1=exp(a_i)*0.5*(t-min(d$Week)) * (as.numeric(wk%*%exp(d_s$v)))
      survi=exp(-z1) 
      
    }
    
    return(survi)
  }
  
  #Bisection for Time to Event.
  bisect=function(i,x,interv)
  {
    mid=as.integer(mean(interv))
    inter=NULL
    sur=Surv_id(mid,i)
    if(x>=sur)
    {
      inter=c(min(interv),mid)
    }
    else
    {
      inter=c(mid,max(interv))
    }
    return(inter)
  }
  
  #Use Bisection method to find out time to event
  Time_to=function(i)
  {
    #Subset Data set
    d=subset(dataL, dataL$id==i)
    #Initialize
    i0=c(as.integer(min(d$Week)*max(Censor)),as.integer(d$Censor[1]*max(Censor)))
    #Bisect Intervals
    while(abs(diff(i0))>=2)
    {
      i0=bisect(i,U[i],i0)
      
    }
    toe=NULL
    #Time to event
    toe=max(t(sapply(min(i0):max(i0),function(x) c(x,Surv_id(x,i))))[,1])
    #If time-to-event ==Censoring time , then take Censor-1
    toe1=ifelse(toe==as.integer(d$Censor[1]*max(Censor)), toe-1,toe)
    return(toe1)
    
  }
  
  #Generate Uniforms
  U=runif(N,0,1)
  
  #Relapsed Indicator
  Surv_Censor=sapply(1:N, function(i) Surv_id(Censor[i],i))
  #if U[i]<S(C[i]) => T[i]>C[i] => Censored => Relapsed =No, T[i]=True(Unobserved) time of relapse
  Relapsed=ifelse(U<Surv_Censor,"No","Yes")
  
  #Relapsed : Yes:1 , No:0
  
  V=sapply(1:N, function(x) if(Relapsed[x]=='Yes') c(x,Time_to(x),1) else c(x,Censor[x],0) )
  
  V1=data.frame(t(V))
  colnames(V1)=c("id",'Time_to_Event','Relapsed')
  V1$Time_to_Event=as.numeric(V1$Time_to_Event)
  V1$id=as.integer(as.numeric(V1$id))
  V1$Censor_week=Censor[V1$id]
  dataL=merge(x=dataL,y=V1, by=c("id"), all.x=T)
  return(dataL)
}

#Creating data list and Init list for QRJM
data_init_list=function(DataSim,Thetas, tau)
{
  #Computation Data list
  max_censor=max(DataSim$Censor_week)
  DataSim$Time_to_Event=DataSim$Time_to_Event/max_censor
  
  DataSim=DataSim[with(DataSim, order(id, Week)), ]
  
  id <- as.integer(transform(DataSim , id=as.numeric(factor(id)))$id)
  offset <- as.vector(c(1, 1+cumsum(tapply(id,id,length))))
  timeVar <- "Week"
  times <-as.vector(DataSim[[timeVar]])
  
  
  #Number of Ids
  n<-length(unique(DataSim$id))
  
  #Longitudinal Betas
  X=DataSim[,c('id',colnames(Thetas$Beta))]
  
  #response
  y=DataSim[,c("id","Week","ANC","PLT")]
  
  #baseline covariates
  W=DataSim[,c("id",names(Thetas$gamma))]
  W=W[!duplicated(W$id),]
  W=W[with(W, order(id)), ]
  rownames(W)=NULL
  
  #Survival Response
  y_Surv <-DataSim%>%
    group_by(id)%>%
    summarise(Time=Time_to_Event[1],
              Event=Relapsed[1])
  y_Surv=y_Surv[with(y_Surv, order(id)), ]
  rownames(y_Surv)=NULL
  
  zeros <- rep(0,nrow(W))

  C <- 10^7
  
  #G-K points
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
  K_pt <- length(sk)   #Points in G-K quadrature
  Time<-y_Surv$Time
  event<-y_Surv$Event
  P <-0.5*(Time-DataSim$Week[offset[c(-length(offset))]])
  
  
  #medicine covariate
  medi=function(id, t)
  {
    d=subset(DataSim, DataSim$id==id)
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
  
  XT=W
  XT$Intercept=1
  XT$Week=Time
  XT$Med1=NA
  XT$Med2=NA
  XT[,c("Med1","Med2")]= t(sapply(1:nrow(XT), function(i) medi(id=XT$id[i],t=XT$Week[i])))
  XT=XT[,colnames(X)]
  
  
  Xs=data.frame(id=rep(XT$id, each=15),Intercept=1,Week=NA,
                Med1=NA,
                Med2=NA,
                sk=rep((sk+1)/2, nrow(XT)),
                stringsAsFactors = FALSE )
  Xs=merge(x=Xs, y=W, by=c("id"), all.x=T )
  Xs=Xs[with(Xs, order(id, Week)), ]
  Xs$toe=rep(Time, each=15)
  Xs$min_week=min(X$Week)
  Xs$Week=(Xs$sk)*(Xs$toe-Xs$min_week)+Xs$min_week
  Xs[,c("Med1","Med2")]= t(sapply(1:nrow(Xs), function(i) medi(id=Xs$id[i],t=Xs$Week[i])))
  Xs=Xs[,colnames(X)]
  
  #Process offset
  offset_Pr=(0:nrow(W))*(length(wk)+1)+1
  
  
  #Longitudinal indices according to Process grid points
  ProL=NULL
  ProS=NULL
  ids=sort(unique(DataSim$id))
  for(i in 1:length(ids))
  {
    d=subset(DataSim, DataSim$id==ids[i])
    tpts=sort(as.vector(c(min(d$Week),(d$Time_to_Event[1]-min(d$Week))*(sk+1)/2 + min(d$Week),d$Time_to_Event[1])))
    d=d[with(d, order(Week)),]
    
    d$index=(i-1)*(length(sk)+1)+sapply(1:nrow(d),  function(j) ifelse(d$Week[j]>=max(tpts),(length(tpts)-1),sum(sapply(2:length(tpts), function(i) ifelse(d$Week[j]>=tpts[i-1]&d$Week[j]<tpts[i], (i-1),0)))))
    tpts1=sqrt(c(min(tpts),diff(tpts[1:(length(tpts)-1)])))
    
    ProL=rbind(ProL,d[,c("id","index")])
    ProS =rbind(ProS,cbind(rep(d$id[1],length(tpts1)),c(1:length(tpts1))+(i-1)*(length(sk)+1),tpts1))
    colnames(ProS)=c("id","Process_time_id","Process_time")
    i=i+1
  }
  ProL=as.matrix(ProL)
  ProS=as.matrix(ProS)
  
  
  initsList_qrjm = list(D_inv = solve(Thetas$Sigma) ,
                        gam.h0 = log(Thetas$lambda_0),
                        lambda=1/Thetas$err_sigma,
                        b = Thetas$Beta ,
                        psi=Thetas$psi ,
                        theta = Thetas$gamma)
  
  #Data list for QRJM
  dat_list_qrjm =list('n'=n,'y'=y,'zeros'=zeros ,'C'=10^7,
                      'K'=length(wk) , 'P'=P , 'wk'=wk , 'offset'=offset,'offset_Pr'=offset_Pr ,'ProL'=ProL,'ProS'=ProS,
                      'X'=X , 'XT'=XT , 'Xs'=Xs , 'ncX'=ncol(X) ,'ncXT'=ncol(XT) ,'ncXs'=ncol(Xs) , 'ncy'=(ncol(y)-2) ,
                      'W'=W , 'ncW'=ncol(W) ,
                      'event'=event ,'ze'=rep(0,(ncol(y)-2)),'Omega' = diag(rep(1,(ncol(y)-2))) ,
                      't1'=as.vector(sapply(1:length(tau), function(t) (1-2*tau[t])/(tau[t]*(1-tau[t])))), #ALD decomposition constant 1
                      't2'=as.vector(sapply(1:length(tau), function(t) 2/(tau[t]*(1-tau[t]))))) #ALD decomposition constant 2) 
 
  
  Grand_list=list('init'=initsList_qrjm, 'data_list'=dat_list_qrjm) 
  return(Grand_list)

}



#-------------------------------------------------------------------------------------------#
#                                       QRJM Model
#-------------------------------------------------------------------------------------------#
sink(paste(path_out,"QRJM.model.txt", sep="/"))
cat("
model {

for(i in 1:n) 
{
 #BM Process
 #----------
 u[offset_Pr[i],1:ncy] ~ dmnorm(ze[],D_inv[,])
 Process[offset_Pr[i],1:ncy] = u[offset_Pr[i],1:ncy] * ProS[offset_Pr[i],3]
 
 for( w in (offset_Pr[i]+1):(offset_Pr[i+1]-1))
 {
 u[w,1:ncy] ~ dmnorm(ze[],D_inv[,])
 Process[w,1:ncy] = Process[(w-1),1:ncy] + u[w,1:ncy] * ProS[w,3]
 }
 
 
 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
 
 for(k in 1:ncy)
    {
    #Exponential variable
    Ex[j,k] ~ dexp(lambda[k])
    
    #mean
    mu[j,k] =  inprod(b[k,1:(ncX-1)], X[j,2:ncX]) +  Process[ProL[j,2],k]  + t1[k]*Ex[j,k]
    tau[j,k] = lambda[k]/(t2[k]*Ex[j,k])
    y[j,(2+k)] ~ dnorm(mu[j,k], tau[j,k])

    } #loop of k
  } #loop of j

 
 
 #Survival part 1
 #---------------
 
 #log.h0T[i] = gam.h0 #Baseline log hazard
 theta.W[i] = inprod(theta[1:(ncW-1)],W[i,2:ncW]) #Baseline Covariate 
 for(k in 1:ncy)
 {
 muT[i,k] = inprod(b[k,1:(ncX-1)],XT[i,2:ncXT]) +  Process[(offset_Pr[i+1]-1),k]  #Fixed part of Longitudinal
 } 
 
 
 #Survival part 2
 #---------------
 
 for(l in 1:K)
  {
  #log.h0s[i,l] = gam.h0 #Baseline log hazard
  for(k in 1:ncy)
 {
 mus[i,k,l] = inprod(b[k,1:(ncX-1)],Xs[(K*(i-1)+l),2:ncXs]) +  Process[((i-1)*(K+1) +1+l),k]  #Fixed part of Longitudinal
 } 
  SurvLong[i,l] = (wk[l]* exp(gam.h0  + inprod(psi[1:ncy], mus[i,1:ncy,l])))

  }

  
  #Zeros trick in JAGS for JM
  #---------------------------------
  phi[i] = C- ((event[i]*(gam.h0 + theta.W[i] + inprod(psi[1:ncy],muT[i,1:ncy]))) -exp(theta.W[i]) * P[i]* sum(SurvLong[i,1:K]))
  zeros[i] ~ dpois(phi[i])
}#loop of i

#Prior
#-------------------------------------
  #Exponential parameter
  #-----------------------------------
  for(k in 1:ncy)
  {
  lambda[k] ~ dgamma(0.001,0.001)
  sig[k] = 1/lambda[k]
  }
  
  
  #Cov matrix for BM Processes
  #-----------------------------------
 
  D_inv[1:ncy,1:ncy] ~ dwish(Omega[,], (ncy+1))
  D[1:ncy,1:ncy] <-inverse(D_inv[1:ncy,1:ncy])
  
  #Beta
  #-----
  #Polynomial and medicine and fixed variables
  #-------------------------------------------
  for(k in 1:ncy)
  {
  for (j in 1:(ncX-1)) {
   b[k,j] ~ dnorm(0,0.001)
  }
  }
  

  #psi
  #-----
  for(k in 1:ncy)
  {
  psi[k] ~ dnorm(0,0.001)
  }

  #theta
  #-----
  for(h in 1:(ncW-1))
  {
  theta[h] ~  dnorm(0,0.001)
  }
  
  
  
  #Prior for penalized B-spline basis coefficients
  #------------------------------------------------
  gam.h0 ~  dmnorm(0,1000)

  
}
",fill = TRUE)
sink()




th1=TP_thetas(N=75, dist.string_ANC = 'ALD0.25', dist.string_PLT = 'ALD0.25')
th2=TP_thetas(N=75, dist.string_ANC = 'ALD0.25', dist.string_PLT = 'ALD0.75')
th3=TP_thetas(N=75, dist.string_ANC = 'ALD0.5', dist.string_PLT = 'ALD0.5')
th4=TP_thetas(N=75, dist.string_ANC = 'ALD0.75', dist.string_PLT = 'ALD0.25')
th5=TP_thetas(N=75, dist.string_ANC = 'ALD0.75', dist.string_PLT = 'ALD0.75')

ds1=Data_Gen(Thetas=th1)
ds1$type="1"
ds2=Data_Gen(Thetas=th2)
ds2$type="2"
ds3=Data_Gen(Thetas=th3)
ds3$type="3"
ds4=Data_Gen(Thetas=th4)
ds4$type="4"
ds5=Data_Gen(Thetas=th5)
ds5$type="5"

ds2$id=ds2$id+length(unique(ds1$Child.Id))
ds3$id=ds3$id+max(ds2$id)
ds4$id=ds4$id+max(ds3$id)
ds5$id=ds5$id+max(ds4$id)
ds=data.frame(rbind(ds1,ds2,ds3,ds4,ds5))
ds=ds[ds$Time_to_Event>30,]
ds$Child.Id=sapply(1:nrow(ds), function(x) paste0('UPN_',ds$id[x],sep=""))


df_names=ds%>%
  group_by(Child.Id)%>%
  summarise(id_x=id[1])

df_names=df_names[with(df_names, order(id_x)),]
df_names$id_new=1:nrow(df_names)
df_names$ChildId_new=sapply(1:nrow(df_names), function(x) paste0('UPN_',df_names$id_new[x],sep=""))
ds=merge(x=ds, y=df_names[,c('Child.Id','id_new','ChildId_new')],by=c('Child.Id'), all.x=T)
ds$Child.Id=ds$ChildId_new
ds$id=ds$id_new
ds=ds[,!colnames(ds)%in%c('id_new','ChildId_new')]
ds=ds[with(ds, order(id)),]



#Checks
u=NULL
u=ds%>%
  group_by(id)%>%
  summarise(Relapsed=Relapsed[1], Censor=Censor_week[1],TOE=Time_to_Event[1],type=type[1])
round(table(u$type)*100/nrow(u),1)
mean(u$Relapsed)

ds=ds[,!colnames(ds)%in%c("type")]


saveRDS(ds,paste(path_out,"/Data/The_Data",".rds",sep="")) #The data
saveRDS(th1,paste(path_out,"/Thetas/Theta1",".rds",sep=""))#(0.25,0.25)
saveRDS(th2,paste(path_out,"/Thetas/Theta2",".rds",sep=""))#(0.25,0.75)
saveRDS(th3,paste(path_out,"/Thetas/Theta3",".rds",sep=""))#(0.5,0.5)
saveRDS(th4,paste(path_out,"/Thetas/Theta4",".rds",sep=""))#(0.75,0.25)
saveRDS(th5,paste(path_out,"/Thetas/Theta5",".rds",sep=""))#(0.75,0.75)



#Set tau
tau=c(0.75,0.75)
Get_names=function(tau)
{
  ths=NULL
  name_qrjm=NULL
  if(tau[1]==0.25 & tau[2]==0.25)
  {
    ths='Theta1'
  }
  if(tau[1]==0.25 & tau[2]==0.75)
  {
    ths='Theta2'
  }
  if(tau[1]==0.5 & tau[2]==0.5)
  {
    ths='Theta3'
  }
  if(tau[1]==0.75 & tau[2]==0.25)
  {
    ths='Theta4'
  }
  if(tau[1]==0.75 & tau[2]==0.75)
  {
    ths='Theta5'
  }
  z=NULL
  name_qrjm=paste('QRJM',as.character(as.integer(tau[1]*100)),as.character(as.integer(tau[2]*100)),sep="_")
  z=c(ths, name_qrjm)
  return(z)
  
}


dsim=readRDS(paste(path_out,"/Data/The_Data",".rds",sep=""))
tn=Get_names(tau)
theta_init=readRDS(paste(path_out,"/Thetas/",tn[1],".rds",sep=""))

#Get the lists
d_init=NULL
d_init=data_init_list(DataSim=dsim,Thetas=theta_init, tau=tau)
#inits1=function(){d_init$init}
inits1=NULL
dat_list=d_init$data_list
Sys.time()
#Run QRJM
jagsfit_qrjm=NULL
jagsfit_qrjm <- jags(data=dat_list, 
                     inits=inits1, 
                     parameters.to.save=c('u','b','D','sig','psi','theta','gam.h0', 'Process'),
                     model.file=paste(path_out,'QRJM.model.txt', sep="/"),parallel = TRUE,
                     n.chains=2, n.iter=11000, n.burnin=1000,
                     n.thin=10,n.adapt = 1000,
                     DIC=TRUE)

saveRDS(jagsfit_qrjm,paste(path_out,"/Results/",tn[2],".rds",sep=""))


