##########################################################################################
## ---------------------------------- Parameter space ------------------------------------
##########################################################################################

cube<-matrix(c(0,-16,-0.5,
               0.5,-16,-0.5,
               0.5,81,-0.5,
               0,81,-0.5,
               0,-16,1,
               0.5,-16,1,
               0.5,81,1,
               0,81,1),ncol=3,byrow=TRUE)
cube<-data.frame(tau=cube[,1],mu=cube[,2],xi=cube[,3])

fig_0 <- plot_ly(type = 'mesh3d',
                 x = cube$tau,
                 y = cube$mu,
                 z = cube$xi,
                 i = c(0,0,3,3,1,1,0,0,0,0,4,4),
                 j = c(1,5,2,6,2,6,3,7,1,2,5,6),
                 k = c(5,4,6,7,6,5,7,4,2,3,6,7),
                 facecolor = c(rep("#48a0cf",4), rep("rgba(221,81,58,1)",4),
                               rep("#cacf48",4))
)%>% layout(title = "Parameter space in 3D",
            scene = list(
              camera = list(eye = list(x = 2.4/1.6, y = -2.7/1.6, z = 1/1.6)),
              xaxis = list(ticktext = list("0", "", "", "∞"), 
                           tickvals = list(0, 0.17, 0.34, 0.5),tickmode = "array",
                           titlefont = list(color="#48a0cf",size = 23),
                           tickfont = list(color="#48a0cf",size = 13),title = "𝜏"),
              yaxis = list(ticktext = list("-∞", "", "", "∞"), 
                           tickvals = list(-16, 17.5, 51, 81),tickmode = "array",
                           titlefont = list(color="rgba(221,81,58,1)",size = 23),
                           tickfont = list(color="rgba(221,81,58,1)",size = 13),title = "𝜇"),
              zaxis = list(ticktext = list("-0.5", "", "", "","","","∞"), 
                           tickvals = list(-0.5, -0.25, 0, 0.25,0.5,0.75,1),tickmode = "array",
                           titlefont = list(color="#cacf48",size = 23),
                           tickfont = list(color="#cacf48",size = 13),title = "𝜉")
            ))%>% style(hoverinfo = 'none')

fig_0




##########################################################################################
## ---------------------------------- 3d support plot-------------------------------------
##########################################################################################

## Positive shape
beta0<-20-0.5/0.2
tau<-seq(0,0.5,length.out = 100)
mu<- seq(81,18,length.out = 100)
tmp<-expand.grid(mu=mu,tau=tau)
tmp<-tmp[tmp$tau==0.5|tmp$mu==81,]
tmp<-rbind(c(18,0),tmp)
tmp<-cbind(tmp,xi=tmp$tau/(tmp$mu-beta0))
num0<-nrow(tmp)-1

xi_t <- 1; tau_t<-seq(0.5,0,length.out = 50); tmp_t<-cbind(mu=beta0+tau_t/xi_t,tau=tau_t,xi=xi_t)
tmp<-rbind(tmp, tmp_t)
num<-nrow(tmp)-1

tmp<-rbind(tmp,c(2*beta0-51, 0, 1))
tmp<-rbind(tmp,c(2*beta0-51, 0.5, 1))
tmp<-rbind(tmp,c(2*beta0-51, 0, 0))
tmp<-rbind(tmp,c(2*beta0-51, 0.5, 0))

## Negative shape
FIRST_STEP<-nrow(tmp)
yn <- 51
mu_N<- seq(-16,yn-1,length.out = 100)
tmp_N<-expand.grid(mu=mu_N,tau=tau)
tmp_N<-tmp_N[tmp_N$tau==0.5|tmp_N$mu==-16,]
tmp_N<-rbind(c(yn,0),tmp_N)
tmp_N<-cbind(tmp_N,xi=tmp_N$tau/(tmp_N$mu-yn))
tmp_N$xi[1]<-0
num0_N<-FIRST_STEP+nrow(tmp_N)-1

xi_t <- -0.5; tau_t<-seq(0.5,0,length.out = 50); tmp_t<-cbind(mu=yn+tau_t/xi_t,tau=tau_t,xi=xi_t)
tmp_N<-rbind(tmp_N, tmp_t)
num_N<-FIRST_STEP+nrow(tmp_N)-1

tmp_N<-rbind(tmp_N,c(81, 0, -0.5))
tmp_N<-rbind(tmp_N,c(81, 0.5, -0.5))
tmp_N<-rbind(tmp_N,c(81, 0, 0))
tmp_N<-rbind(tmp_N,c(81, 0.5, 0))
tmp<-rbind(tmp,tmp_N)



fig <- plot_ly(type = 'mesh3d',
               x = tmp$tau,
               y = tmp$mu,
               z = tmp$xi,
               i = c(rep(0,num-1),num0,num0,num+3,num+3,num+3,num+4,rep(num+4,100),num+3,num+3,
                     rep(FIRST_STEP,num-1), num0_N, num0_N, num_N+3,num_N+3,num_N+3,num_N+4,rep(num_N+4,100),num_N+4,num_N+4,
                     num_N+3,num_N+3),
               j = c(1:(num-1), num, num+1, num+1,num+2,num+4,num+2,num0:(num0-99),num+1,num,
                     (FIRST_STEP+1):(num_N-1),num_N,num_N+1,num_N+1,num_N+2,num_N+4,num_N+2,num0_N:(num0_N-99),FIRST_STEP+100,num+4,
                     num_N+1,num_N),
               k = c(2:num,num+1,num+2,     num+2,num+4,FIRST_STEP+100,num0,(num0-1):(num0-100),num,0,
                     (FIRST_STEP+2):num_N, num_N+1, num_N+2,num_N+2,num_N+4,100,    num0_N,(num0_N-1):(num0_N-100),num+4,100,
                     num_N,FIRST_STEP),
               #intensity = c(0, 0.33, 0.66, 1),
               facecolor = c(rep("#1a1c1c",num-1), rep("#cacf48",2),
                             rep("#48a0cf",3),rep("rgba(221,81,58,1)",101),
                             rep("rgba(221,81,58,1)",2),
                             rep("#1a1c1c",num-1), rep("#cacf48",2),
                             rep("#48a0cf",3),rep("rgba(221,81,58,1)",103),
                             rep("rgba(221,81,58,1)",2))
               #colors = colorRamp(c("red", "green", "blue"))
)%>% layout(title = "GEV support in 3D",
            scene = list(
              camera = list(eye = list(x = 2.4/1.6, y = -2.7/1.6, z = 1/1.6)),
              xaxis = list(ticktext = list("0", "", "", "∞"), 
                           tickvals = list(0, 0.17, 0.34, 0.5),tickmode = "array",
                           titlefont = list(color="#48a0cf",size = 23),
                           tickfont = list(color="#48a0cf",size = 13),title = "𝜏"),
              yaxis = list(ticktext = list("-∞", "Y(1)", "Y(n)", "∞"), 
                           tickvals = list(-16, 17.5, 51, 81),tickmode = "array",
                           titlefont = list(color="rgba(221,81,58,1)",size = 23),
                           tickfont = list(color="rgba(221,81,58,1)",size = 13),title = "𝜇"),
              zaxis = list(ticktext = list("-0.5", "", "0", "","","","∞"), 
                           tickvals = list(-0.5, -0.25, 0, 0.25,0.5,0.75,1),tickmode = "array",
                           titlefont = list(color="#cacf48",size = 23),
                           tickfont = list(color="#cacf48",size = 13),title = "𝜉")
            ))%>% style(hoverinfo = 'none')

fig
Sys.setenv("plotly_username"="lfz5044")
Sys.setenv("plotly_api_key"="nn8hSl9HkAYX9RMvHw3L")
chart_link = api_create(fig, filename="support")
chart_link




##########################################################################################
## ---------------------------------- Cross sections -------------------------------------
##########################################################################################

## Positive shape
h<-0.1
x0<-c(0,0.5);y0<-c(-16,81);z0<-matrix(0,ncol=2,nrow=2);z0_1<-matrix(h,ncol=2,nrow=2)
xl<-c(0,0.5);yl<-c(beta0,beta0+0.5/h);zl<-c(h,h)

fig <- plot_ly()%>%add_trace(opacity = 0.91,type='mesh3d',
               x = tmp$tau,
               y = tmp$mu,
               z = tmp$xi,
               i = c(rep(0,num-1),num0,num0,num+3,num+3,num+3,num+4,rep(num+4,100),num+3,num+3,
                     rep(FIRST_STEP,num-1), num0_N, num0_N, num_N+3,num_N+3,num_N+3,num_N+4,rep(num_N+4,100),num_N+4,num_N+4,
                     num_N+3,num_N+3),
               j = c(1:(num-1), num, num+1, num+1,num+2,num+4,num+2,num0:(num0-99),num+1,num,
                     (FIRST_STEP+1):(num_N-1),num_N,num_N+1,num_N+1,num_N+2,num_N+4,num_N+2,num0_N:(num0_N-99),FIRST_STEP+100,num+4,
                     num_N+1,num_N),
               k = c(2:num,num+1,num+2,     num+2,num+4,FIRST_STEP+100,num0,(num0-1):(num0-100),num,0,
                     (FIRST_STEP+2):num_N, num_N+1, num_N+2,num_N+2,num_N+4,100,    num0_N,(num0_N-1):(num0_N-100),num+4,100,
                     num_N,FIRST_STEP),
               facecolor = c(rep("rgba(86, 87, 75,0.5)",num-1), rep("rgba(86, 87, 75,0.5)",2),
                             rep("rgba(86, 87, 75,0.5)",3),rep("rgba(86, 87, 75,0.5)",101),
                             rep("rgba(86, 87, 75,0.5)",2),
                             rep("rgba(86, 87, 75,0.5)",num-1), rep("rgba(86, 87, 75,0.5)",2),
                             rep("rgba(86, 87, 75,0.5)",3),rep("rgba(86, 87, 75,0.5)",103),
                             rep("rgba(86, 87, 75,0.5)",2))
   )%>% add_trace(type="surface",z=~z0_1,x=~x0,y=~y0,surfacecolor=z0_1,colorscale=list(c(0,1),c("steelblue","steelblue")),showscale=FALSE)%>% 
        add_trace(x=~xl,y=~yl,z=~zl,type='scatter3d',line=list(color = 'black', width = 6),mode='lines',show)%>%
        layout(title = "GEV support in 3D",showlegend = FALSE,
            scene = list(
              camera = list(eye = list(x = 2.4/1.6, y = -2.7/1.6, z = 1/1.6)),
              xaxis = list(ticktext = list("0", "", "", "∞"), 
                           tickvals = list(0, 0.17, 0.34, 0.5),tickmode = "array",
                           titlefont = list(color="#48a0cf",size = 23),
                           tickfont = list(color="#48a0cf",size = 13),range = c(0,0.5),title = "𝜏"),
              yaxis = list(ticktext = list("-∞", "Y(1)", "Y(n)", "∞"), 
                           tickvals = list(-16, 17.5, 51, 81),tickmode = "array",
                           titlefont = list(color="rgba(221,81,58,1)",size = 23),
                           tickfont = list(color="rgba(221,81,58,1)",size = 13),range = c(-16,81),title = "𝜇"),
              zaxis = list(ticktext = list("-0.5", "", "0", "","","","∞"), 
                           tickvals = list(-0.5, -0.25, 0, 0.25,0.5,0.75,1),tickmode = "array",
                           titlefont = list(color="#cacf48",size = 23),
                           tickfont = list(color="#cacf48",size = 13),range = c(-0.5,1),title = "𝜉")
            ))%>% style(hoverinfo = 'none')

fig

## Negative shape
h<- -0.15
x0<-c(0,0.5);y0<-c(-16,81);z0<-matrix(0,ncol=2,nrow=2);z0_1<-matrix(h,ncol=2,nrow=2)
xl<-c(0,0.5);yl<-c(yn,yn+0.5/h);zl<-c(h,h)

fig <- plot_ly()%>%add_trace(opacity = 0.91,type='mesh3d',
                             x = tmp$tau,
                             y = tmp$mu,
                             z = tmp$xi,
                             i = c(rep(0,num-1),num0,num0,num+3,num+3,num+3,num+4,rep(num+4,100),num+3,num+3,
                                   rep(FIRST_STEP,num-1), num0_N, num0_N, num_N+3,num_N+3,num_N+3,num_N+4,rep(num_N+4,100),num_N+4,num_N+4,
                                   num_N+3,num_N+3),
                             j = c(1:(num-1), num, num+1, num+1,num+2,num+4,num+2,num0:(num0-99),num+1,num,
                                   (FIRST_STEP+1):(num_N-1),num_N,num_N+1,num_N+1,num_N+2,num_N+4,num_N+2,num0_N:(num0_N-99),FIRST_STEP+100,num+4,
                                   num_N+1,num_N),
                             k = c(2:num,num+1,num+2,     num+2,num+4,FIRST_STEP+100,num0,(num0-1):(num0-100),num,0,
                                   (FIRST_STEP+2):num_N, num_N+1, num_N+2,num_N+2,num_N+4,100,    num0_N,(num0_N-1):(num0_N-100),num+4,100,
                                   num_N,FIRST_STEP),
                             facecolor = c(rep("rgba(86, 87, 75,0.5)",num-1), rep("rgba(86, 87, 75,0.5)",2),
                                           rep("rgba(86, 87, 75,0.5)",3),rep("rgba(86, 87, 75,0.5)",101),
                                           rep("rgba(86, 87, 75,0.5)",2),
                                           rep("rgba(86, 87, 75,0.5)",num-1), rep("rgba(86, 87, 75,0.5)",2),
                                           rep("rgba(86, 87, 75,0.5)",3),rep("rgba(86, 87, 75,0.5)",103),
                                           rep("rgba(86, 87, 75,0.5)",2))
)%>% add_trace(type="surface",z=~z0_1,x=~x0,y=~y0,surfacecolor=z0_1,colorscale=list(c(0,1),c("steelblue","steelblue")),showscale=FALSE)%>% 
  add_trace(x=~xl,y=~yl,z=~zl,type='scatter3d',line=list(color = 'black', width = 6),mode='lines')%>%
  layout(title = "GEV support in 3D",showlegend = FALSE,
         scene = list(
           camera = list(eye = list(x = 2.4/1.6, y = -2.7/1.6, z = 1/1.6)),
           xaxis = list(ticktext = list("0", "", "", "∞"), 
                        tickvals = list(0, 0.17, 0.34, 0.5),tickmode = "array",
                        titlefont = list(color="#48a0cf",size = 23),
                        tickfont = list(color="#48a0cf",size = 13),range = c(0,0.5),title = "𝜏"),
           yaxis = list(ticktext = list("-∞", "Y(1)", "Y(n)", "∞"), 
                        tickvals = list(-16, 17.5, 51, 81),tickmode = "array",
                        titlefont = list(color="rgba(221,81,58,1)",size = 23),
                        tickfont = list(color="rgba(221,81,58,1)",size = 13),range = c(-16,81),title = "𝜇"),
           zaxis = list(ticktext = list("-0.5", "", "0", "","","","∞"), 
                        tickvals = list(-0.5, -0.25, 0, 0.25,0.5,0.75,1),tickmode = "array",
                        titlefont = list(color="#cacf48",size = 23),
                        tickfont = list(color="#cacf48",size = 13),range = c(-0.5,1),title = "𝜉")
         ))%>% style(hoverinfo = 'none')

fig



##########################################################################################
## ----------------------------------- Contour plots -------------------------------------
##########################################################################################

# calculate the likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  if(xi==0) return(-length(Y)*log(tau)-sum((Y-mu)/tau)-sum(exp(-(Y-mu)/tau)))
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)
beta_hat_line <- function(beta_hat,tau, xi) beta_hat + tau/xi

set.seed(123)
xi_try<-0.2
mu <- 20
tau <- 0.5
n <- 1000
Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
delta<-0.02
mu+delta-(tau+delta)/(xi_try-delta)
Res <- optim(c(mu+delta,xi_try+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi_try-delta,tau-delta), upper = c(mu+delta,xi_try+delta,tau+delta), control=list(fnscale=-1))
beta_hat <- Res$par[1]-Res$par[3]/Res$par[2]


m<-1.5  # Controls the cross section: xi_cross = m*xi_MLE
resolution <- 250
xi_MLE <- Res$par[2]
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(16,24.5,length.out = resolution)
Par=expand.grid(tau=Tau,mu=Mu)
Par<-cbind(mu=Par$mu, xi=m*xi_MLE, tau=Par$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)
volcano <- t(matrix(lik_Par,length(Tau)))


beta_t <- seq(14.4,18.5,length.out = 12)
xi_t<-m*xi_MLE
#---1
lik_Beta1<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[1]+Tau[i]/xi_t
  lik_Beta1[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta1)<-c('tau','mu','Log_likelihood')
lik_Beta1<-data.frame(lik_Beta1[-1,])
lik_Beta1[lik_Beta1$Log_likelihood< -12, 3]<- -12

#---2
lik_Beta2<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[2]+Tau[i]/xi_t
  lik_Beta2[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta2)<-c('tau','mu','Log_likelihood')
lik_Beta2<-data.frame(lik_Beta2[-1,])
lik_Beta2[lik_Beta2$Log_likelihood< -12, 3]<- -12

#---3
lik_Beta3<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[3]+Tau[i]/xi_t
  lik_Beta3[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta3)<-c('tau','mu','Log_likelihood')
lik_Beta3<-data.frame(lik_Beta3[-1,])
lik_Beta3[lik_Beta3$Log_likelihood< -12, 3]<- -12

#---4
lik_Beta4<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[4]+Tau[i]/xi_t
  lik_Beta4[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta4)<-c('tau','mu','Log_likelihood')
lik_Beta4<-data.frame(lik_Beta4[-1,])
lik_Beta4[lik_Beta4$Log_likelihood< -12, 3]<- -12

#---5
lik_Beta5<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[5]+Tau[i]/xi_t
  lik_Beta5[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta5)<-c('tau','mu','Log_likelihood')
lik_Beta5<-data.frame(lik_Beta5[-1,])
lik_Beta5[lik_Beta5$Log_likelihood< -12, 3]<- -12

#---6
lik_Beta6<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[6]+Tau[i]/xi_t
  lik_Beta6[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta6)<-c('tau','mu','Log_likelihood')
lik_Beta6<-data.frame(lik_Beta6[-1,])
lik_Beta6[lik_Beta6$Log_likelihood< -12, 3]<- -12

#---7
lik_Beta7<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[7]+Tau[i]/xi_t
  lik_Beta7[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta7)<-c('tau','mu','Log_likelihood')
lik_Beta7<-data.frame(lik_Beta7[-1,])
lik_Beta7[lik_Beta7$Log_likelihood< -12, 3]<- -12

#---8
lik_Beta8<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[8]+Tau[i]/xi_t
  lik_Beta8[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta8)<-c('tau','mu','Log_likelihood')
lik_Beta8<-data.frame(lik_Beta8[-1,])
lik_Beta8[lik_Beta8$Log_likelihood< -12, 3]<- -12

#---9
lik_Beta9<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[9]+Tau[i]/xi_t
  lik_Beta9[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta9)<-c('tau','mu','Log_likelihood')
lik_Beta9<-data.frame(lik_Beta9[-1,])
lik_Beta9[lik_Beta9$Log_likelihood< -12, 3]<- -12

#---10
lik_Beta10<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[10]+Tau[i]/xi_t
  lik_Beta10[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta10)<-c('tau','mu','Log_likelihood')
lik_Beta10<-data.frame(lik_Beta10[-1,])
lik_Beta10[lik_Beta10$Log_likelihood< -12, 3]<- -12

#---11
lik_Beta11<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[11]+Tau[i]/xi_t
  lik_Beta11[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta11)<-c('tau','mu','Log_likelihood')
lik_Beta11<-data.frame(lik_Beta11[-1,])
lik_Beta11[lik_Beta11$Log_likelihood< -12, 3]<- -12

#---12
lik_Beta12<-matrix(NA,nrow=length(Tau),ncol=3)
for(i in 1:length(Tau)){
  mu_t<- beta_t[12]+Tau[i]/xi_t
  lik_Beta12[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
}
colnames(lik_Beta12)<-c('tau','mu','Log_likelihood')
lik_Beta12<-data.frame(lik_Beta12[-1,])
lik_Beta12[lik_Beta12$Log_likelihood< -12, 3]<- -12

#---profile
beta_T <- sort(c(seq(13,19,length.out = 30),beta_t))
lik_Beta<-matrix(NA,nrow=length(Tau),ncol=3)
prof <- matrix(NA,nrow=length(beta_T),ncol=3)
for(j in 1:length(beta_T)){
  for(i in 1:length(Tau)){
    mu_t<- beta_T[j]+Tau[i]/xi_t
    lik_Beta[i,]<-c(Tau[i],mu_t,Lik(c(mu_t,xi_t,Tau[i]),Y)/1000)
  }
  prof[j,]<-lik_Beta[which.max(lik_Beta[,3]),]
}

colnames(prof)<-c('tau','mu','Log_likelihood')
prof<-data.frame(prof)


z1<-matrix(-10,resolution,resolution)
## Plot the log-lik evaluated
fig <- plot_ly(lik_Beta1, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', 
               mode = 'lines', colors = terrain.colors(40), color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta2, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta3, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta4, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta5, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta6, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta7, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta8, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta9, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta10, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta11, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta12, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(type='surface', x=~Tau, y=~Mu, z=~z1, surfacecolor=volcano)%>%
  hide_colorbar() %>%
  layout(title="Log-likelihood on one cross section", showlegend = FALSE,
         scene = list(
           camera = list(eye = list(x = 2.4/1.5, y = -2.7/1.5, z = 1/1.5)),
           xaxis = list(ticktext = list("0", "", "", "∞"), 
                        tickvals = list(0, 0.5, 1, 1.5),tickmode = "array",
                        titlefont = list(color="#48a0cf",size = 23),
                        tickfont = list(color="#48a0cf",size = 13),range = c(0,1.5),title = "𝜏"),
           yaxis = list(ticktext = list("-∞", "Y(1)", "","", "∞"), 
                        tickvals = list(16.425, min(Y), 21.942, 24.702, 27.46),tickmode = "array",
                        titlefont = list(color="rgba(221,81,58,1)",size = 23),
                        tickfont = list(color="rgba(221,81,58,1)",size = 13),range = c(16.425,27.46),title = "𝜇"),
           zaxis = list(#ticktext = list("", "","","0"), 
             #tickvals = list(-15, -10, -5 ,0),tickmode = "array", 
             titlefont = list(color="#cacf48",size = 13),
             tickfont = list(color="#cacf48",size = 9),range = c(-10,0),title = "Log-lik")
         ))%>% style(hoverinfo = 'none')

fig


xl<-rep(prof[which.max(prof$Log_likelihood),1],2)
yl<-rep(prof[which.max(prof$Log_likelihood),2],2)
zl<-c(-10,prof[which.max(prof$Log_likelihood),3])

Tmp<-data.frame(xl=xl,yl=yl,Log_likelihood=zl)
fig <- plot_ly(lik_Beta1, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', 
               mode = 'lines', colors = terrain.colors(40), color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta2, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta3, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta4, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta5, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta6, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta7, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta8, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta9, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta10, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta11, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=lik_Beta12, x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',color = ~Log_likelihood)%>%
  add_trace(data=prof,       x = ~tau, y = ~mu, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',line=list(color = 'blue', width = 4))%>%
  add_trace(data=Tmp, x = ~xl, y = ~yl, z = ~Log_likelihood, type = 'scatter3d', mode = 'lines',line=list(color = 'gray', width = 3, dash = 'dash'))%>%
  add_trace(type='surface', x=~Tau, y=~Mu, z=~z1, surfacecolor=volcano)%>%
  hide_colorbar() %>%
  layout(title="Log-likelihood on one cross section", showlegend = FALSE,
         scene = list(
           camera = list(eye = list(x = 2.4/1.5, y = -2.7/1.5, z = 1/1.5)),
           xaxis = list(ticktext = list("0", "", "", "∞"), 
                        tickvals = list(0, 0.5, 1, 1.5),tickmode = "array",
                        titlefont = list(color="#48a0cf",size = 23),
                        tickfont = list(color="#48a0cf",size = 13),range = c(0,1.5),title = "𝜏"),
           yaxis = list(ticktext = list("-∞", "Y(1)", "","", "∞"), 
                        tickvals = list(16.425, min(Y), 21.942, 24.702, 27.46),tickmode = "array",
                        titlefont = list(color="rgba(221,81,58,1)",size = 23),
                        tickfont = list(color="rgba(221,81,58,1)",size = 13),range = c(16.425,27.46),title = "𝜇"),
           zaxis = list(#ticktext = list("", "","","0"), 
                        #tickvals = list(-15, -10, -5 ,0),tickmode = "array", 
                      titlefont = list(color="#cacf48",size = 13),
                        tickfont = list(color="#cacf48",size = 9),range = c(-10,0),title = "Log-lik")
         ))%>% style(hoverinfo = 'none')

fig



##########################################################################################
## ----------------------------- Multiple cross sections ---------------------------------
##########################################################################################

M<-seq(0.2,1.8,by=0.4)  # Controls the cross section: xi_cross = m*xi_MLE
resolution <- 250
xi_MLE <- Res$par[2]
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(16,27.46,length.out = resolution)
Par0=expand.grid(tau=Tau,mu=Mu)
Points<- data.frame(NA,ncol=3,nrow=length(M))
colnames(Points) <- c('mu','xi','tau')

## --1
Par<-cbind(mu=Par0$mu, xi=M[1]*xi_MLE,tau=Par0$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp/1000
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)
volcano1 <- t(matrix(lik_Par,length(Tau)))
z1<-matrix(M[1]*xi_MLE,resolution,resolution)
Points[1,] <- Par[which.max(lik_Par),]
  
  
## --2
Par<-cbind(mu=Par0$mu, xi=M[2]*xi_MLE,tau=Par0$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp/1000
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)
volcano2 <- t(matrix(lik_Par,length(Tau)))
z2<-matrix(M[2]*xi_MLE,resolution,resolution)
Points[2,] <- Par[which.max(lik_Par),]

## --3
Par<-cbind(mu=Par0$mu, xi=M[3]*xi_MLE,tau=Par0$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp/1000
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)
volcano3 <- t(matrix(lik_Par,length(Tau)))
z3<-matrix(M[3]*xi_MLE,resolution,resolution)
Points[3,] <- Par[which.max(lik_Par),]

## --4
Par<-cbind(mu=Par0$mu, xi=M[4]*xi_MLE,tau=Par0$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp/1000
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)
volcano4 <- t(matrix(lik_Par,length(Tau)))
z4<-matrix(M[4]*xi_MLE,resolution,resolution)
Points[4,] <- Par[which.max(lik_Par),]

## --5
Par<-cbind(mu=Par0$mu, xi=M[5]*xi_MLE,tau=Par0$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp/1000
}
range(lik_Par,na.rm = TRUE)
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)
volcano5 <- t(matrix(lik_Par,length(Tau)))
z5<-matrix(M[5]*xi_MLE,resolution,resolution)
Points[5,] <- Par[which.max(lik_Par),]

col_scale<-list(seq(-0,1,length.out=40),terrain.colors(40))
fig <- plot_ly()%>%
  add_surface(opacity = 0.92, x = ~Tau, y = ~Mu, z = ~z1,  cmin = -13, cmax = 0, 
              colorscale = col_scale,surfacecolor = volcano1)%>%
  add_surface(opacity = 0.92, x = ~Tau, y = ~Mu, z = ~z2,  cmin = -11, cmax = 0, 
              colorscale = col_scale,surfacecolor = volcano2)%>%
  add_surface(opacity = 0.92, x = ~Tau, y = ~Mu, z = ~z3,  cmin = -10, cmax = 0, 
              colorscale = col_scale,surfacecolor = volcano3)%>%
  add_surface(opacity = 0.92, x = ~Tau, y = ~Mu, z = ~z4,  cmin = -8, cmax = 0, 
              colorscale = col_scale,surfacecolor = volcano4)%>%
  add_surface(opacity = 0.92, x = ~Tau, y = ~Mu, z = ~z5,  cmin = -7, cmax = 0, 
              colorscale = col_scale,surfacecolor = volcano5)%>%
  add_trace(type='scatter3d',mode='markers', data=Points,x=~tau,y=~mu,z=~xi, marker = list(size = 2,color="darkblue"))%>%
  hide_colorbar() %>%
  layout(title="Log-likelihood on multiple cross sections", showlegend = FALSE,
         scene = list(
           camera = list(eye = list(x = -5.5/2.4, y = -2/2.45, z = 0.89/2.4)),
           xaxis = list(ticktext = list("0", "", "", "∞"), 
                        tickvals = list(0, 0.5, 1, 1.5),tickmode = "array",
                        titlefont = list(color="#48a0cf",size = 23),
                        tickfont = list(color="#48a0cf",size = 13),range = c(0,1.5),title = "𝜏"),
           yaxis = list(ticktext = list("-∞", "Y(1)", "","", "∞"), 
                        tickvals = list(16.425, min(Y), 21.942, 24.702, 27.46),tickmode = "array",
                        titlefont = list(color="rgba(221,81,58,1)",size = 23),
                        tickfont = list(color="rgba(221,81,58,1)",size = 13),range = c(16.425,27.46),title = "𝜇"),
           zaxis = list(#ticktext = list("", "","","0"), 
             #tickvals = list(-15, -10, -5 ,0),tickmode = "array", 
             titlefont = list(color="#cacf48",size = 23),
             tickfont  = list(color="#cacf48",size = 9),range = c(0,0.6),title = "𝜉")
         ))%>% style(hoverinfo = 'none')

fig
# chart_link = api_create(fig, filename="multiple_CS")
# chart_link



##########################################################################################
## ---------------------------------- Convergence rate -----------------------------------
##########################################################################################

# calculate the likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  if(xi==0) return(-length(Y)*log(tau)-sum((Y-mu)/tau)-sum(exp(-(Y-mu)/tau)))
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)
xi_try<-0.2
mu <- 20
tau <- 0.5
resolution <- 100
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(18.5,21.5,length.out = resolution)
xlim=c(-0.3,1.5)
ylim=c(range(Mu)[1]+0.05*diff(range(Mu)),range(Mu)[2]-0.05*diff(range(Mu)))
Par0=expand.grid(tau=Tau,mu=Mu)

N<-c(seq(5,100,by=2),110,150,300,1000,5000,10000,50000)
trace_min<-data.frame(a=rep(NA,length(N)),b=rep(NA,length(N)))
for(iter in 1:length(N)){
  n <- N[iter]
  set.seed(123); Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
  xi0 <- xi_try
  Par<-cbind(mu=Par0$mu, xi=xi0, tau=Par0$tau)
  lik_Par<-rep(NA,nrow(Par))
  for(i in 1:nrow(Par)){
    tmp<-Lik(Par[i,],Y)
    lik_Par[i]<-tmp
  }
  lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)
  
  path <- sprintf("~/Desktop/Research/Thesis/traces/trace%d.png",iter)
  png(path,width=580, height=580)
  plot(x=0.5,y=21,type='n',xlab='',ylab='',xlim=xlim, ylim=ylim)
  title(xlab=expression(tau), ylab=expression(mu), cex.lab=1.5)
  abline(v=0,col="grey",lty=2)
  image(x=Tau,y=Mu,
        z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)
  trace_min[iter,]<-c(min(Y),(min(Y)-bound(-1000,xi0))/1000)
  for(i in 1:iter){
    abline(a=trace_min[i,]$a,b=trace_min[i,]$b,col="grey")
  }
  lines(x=c(0,0),y=c(-100,min(Y)),col='grey')
  points(Par[which.max(lik_Par),][3], Par[which.max(lik_Par),][1], col="blue", pch=20, cex=1.6)
  points(tau,mu, col="red", pch=20, cex=1.9)
  points(tau,mu, col="red",  cex=1.9)
  text(0.02,min(Y),expression(Y[(1)]),pos=2,cex=1.5)
  legend("topleft",cex=1.5,pch=c(20,20), col=c('red','blue'),legend=c(expression(theta[0]),expression(hat(theta)[n])))
  dev.off()
  cat('iter=',iter,'\n')
}


## For negative shape
xi_try<- -0.2
mu <- 20
tau <- 0.5
resolution <- 200
Tau=seq(0,1.7,length.out = resolution)
Mu=seq(18.5,22.5,length.out = resolution)
xlim=c(-0.3,1.5)
ylim=c(19,22)
Par0=expand.grid(tau=Tau,mu=Mu)

n <- 5
set.seed(123); Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
xi0 <- xi_try
Par<-cbind(mu=Par0$mu, xi=xi0, tau=Par0$tau)
lik_Par<-rep(NA,nrow(Par))
for(i in 1:nrow(Par)){
  tmp<-Lik(Par[i,],Y)
  lik_Par[i]<-tmp
}
lik_Par[which(lik_Par<quantile(lik_Par,0.25,na.rm = TRUE))]<-quantile(lik_Par,0.25,na.rm = TRUE)

path <- "~/Desktop/Research/Thesis/traces/trace_neg1.png"
png(path,width=580, height=580)
plot(x=0.5,y=21,type='n',xlab='',ylab='',xlim=xlim, ylim=ylim)
title(xlab=expression(tau), ylab=expression(mu), cex.lab=1.5)
abline(v=0,col="grey",lty=2)
image(x=Tau,y=Mu,
      z=matrix(lik_Par,length(Tau)),col=terrain.colors(40),add=TRUE)
abline(v=0)
lines(x=c(0,0),y=c(-100,max(Y)))
lines(x=c(0,10),y=c(max(Y),max(Y)+10/xi_try))
abline(a=max(Y),b=1/xi_try,col='gray',lty=2)
points(Par[which.max(lik_Par),][3], Par[which.max(lik_Par),][1], col="blue", pch=20, cex=1.6)
points(tau,mu, col="red", pch=20, cex=1.9)
points(tau,mu, col="red",  cex=1.9)
text(0.02,max(Y),expression(Y[(n)]),pos=2,cex=1.5)
legend("topleft",cex=1.5,pch=c(20,20), col=c('red','blue'),legend=c(expression(theta[0]),expression(hat(theta)[n])))
dev.off()



##########################################################################################
## ----------------------------------- Generic plots -------------------------------------
##########################################################################################

# calculate the likelihhood
Lik<-function(par,Y){
  mu<-par[1]; xi<-par[2]; tau<-par[3]
  if(tau<=0) return(NA)
  if(xi==0) return(-length(Y)*log(tau)-sum((Y-mu)/tau)-sum(exp(-(Y-mu)/tau)))
  W <- 1+ (xi)*(Y-mu)/(tau)
  if(any(W<0)) return(NA) else return(-length(Y)*log(tau)-(xi+1)*sum(log(W))/xi-sum(W^{-1/xi}))
}


bound <- function(tau,xi,alpha=0.1) tau/xi+min(Y)
bound_max <- function(tau,xi,alpha=0.1) tau/xi+max(Y)
beta_hat_line <- function(beta_hat,tau, xi) beta_hat + tau/xi

## Positive shape
set.seed(123)
xi_try<-0.2
mu <- 20
tau <- 0.5
n <- 1000
Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
delta<-0.02
mu+delta-(tau+delta)/(xi_try-delta)
Res <- optim(c(mu+delta,xi_try+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi_try-delta,tau-delta), upper = c(mu+delta,xi_try+delta,tau+delta), control=list(fnscale=-1))
beta_hat <- Res$par[1]-Res$par[3]/Res$par[2]

### -- Plot the slice with fixed betas
png("~/Desktop/Research/Thesis/sliced.png",width=485, height=450)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,min(Y),bound(1000,0.1),-1000))
plot(x,type='n',xlab=expression(tau), ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6, cex.lab=1.5, yaxt='n', xaxt='n')
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,min(Y)-0.05*min(Y),expression(Y[(1)]),pos=3,cex=1.5)

alpha=0.7
polygon(x,col=scales::alpha('steelblue',alpha))
z1<-matrix(-10,resolution,resolution)
dev.off()

png("~/Desktop/Research/Thesis/fix_beta.png",width=485, height=450)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(-1000,-1000,min(Y),bound(1000,0.1),-1000))
plot(x,type='n',xlab=expression(tau), ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6, cex.lab=1.5, yaxt='n', xaxt='n')
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,min(Y)-0.05*min(Y),expression(Y[(1)]),pos=3,cex=1.5)

alpha=0.7
polygon(x,col=scales::alpha('steelblue',alpha))
z1<-matrix(-10,resolution,resolution)

betas_can<- seq(min(Y),-10,by=-2)
for(i in 1:length(betas_can)){
  lines(x=c(0,10),y=c(betas_can[i],betas_can[i]+10/0.1))
  lines(x=c(0,-10),y=c(betas_can[i],betas_can[i]-10/0.1),lty=2,col='grey')
}
dev.off()


## Negative shape
set.seed(123)
xi_try<- -0.2
mu <- 20
tau <- 0.5
n <- 1000
Y <- extRemes::revd(n, loc=mu, scale=tau, shape=xi_try)
delta<-0.02
mu+delta-(tau+delta)/(xi_try-delta)
Res <- optim(c(mu+delta,xi_try+delta,tau+delta), Lik, Y=Y, method = "L-BFGS-B", lower = c(mu-delta,xi_try-delta,tau-delta), upper = c(mu+delta,xi_try+delta,tau+delta), control=list(fnscale=-1))
beta_hat <- Res$par[1]-Res$par[3]/Res$par[2]

### -- Plot the slice with fixed betas
png("~/Desktop/Research/Thesis/fix_beta_neg.png",width=485, height=450)
x<-cbind(tau=c(1000,0,0,1000,1000),mu=c(1000,1000,max(Y),bound_max(1000,-0.15),1000))
plot(x,type='n',xlab=expression(tau),ylab='',xlim=c(-0.3,1.5),ylim=c(0.3*min(Y),2.3*min(Y)),
     cex.main=1.6,cex.lab=1.5,yaxt='n', xaxt='n')
title(ylab=expression(mu),line=0, cex.lab=1.5)
axis(side=1, at=c(0,1.565), labels=c('0',expression(infinity)),cex.axis=1.5)
abline(v=0,col="grey",lty=2)
text(-0.12,max(Y)-0.05*max(Y),expression(Y[(n)]),pos=3,cex=1.5)

alpha=0.7
polygon(x,col=scales::alpha('steelblue',alpha))

betas_can<- seq(max(Y),55,by=2)
for(i in 1:length(betas_can)){
  lines(x=c(0,10),y=c(betas_can[i],betas_can[i]-10/0.15))
  lines(x=c(0,-10),y=c(betas_can[i],betas_can[i]+10/0.15),lty=2,col='grey')
}
dev.off()
