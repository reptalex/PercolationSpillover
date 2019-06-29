library(ggplot2)
library(tidyverse)
library(data.table)
library(viridis)
library(ggpubr)
library(mgcv)
library(parallel)
library(segmented)
set.seed(1)

n=1000
rate <- function(a,b,x) exp(a+b*x)

PoisMix <- function(sd.x=1,tot=0,n=1000,mx=6,b=3,d=-3,length.out=50){
  x <- rnorm(n,sd = sd.x)
  a <- seq((-mx)/2,(mx)/2,length.out = length.out)+tot/2
  c <- rev(a)

  D <- data.table(a,b,c,d)
  
  Y <- apply(D,1,FUN=function(m,x) rpois(n = length(x),lambda=rate(m[1],m[2],x)+rate(m[3],m[4],x)),x=x)
  
  fits <- apply(Y,2,FUN=function(y,x) glm(y~x,family=poisson),x=x)
  
  D[,slope:=sapply(fits,coef)[2,]]
  D[,se:=sapply(fits,FUN=function(fit) summary(fit)$coefficients['x','Std. Error'])]
  D[,diff:=a-c]
  D[,total:=tot]
  D[,sd:=sd.x]
  return(D)
}

SDs <- c(0.1,.25,.5,1)

Ds <- lapply(SDs,PoisMix,n=1e3) %>% rbindlist
Ds[,total:=factor(total)]
Ds[,sd:=factor(sd)]
ggplot(Ds,aes(diff,slope,group=sd,color=sd))+
  geom_line(lwd=2)+
  scale_color_manual(values=viridis(length(SDs)))


SDs <- rep(1,500)
D <- lapply(SDs,PoisMix) %>% rbindlist

D[,beta:=6/(1+exp(-diff))-3]
g1 <- ggplot(D,aes(diff,slope))+
  geom_point(cex=5,color=rgb(0,0,0,0.01),pch=16)+
  geom_smooth(lwd=2)+
  geom_line(aes(diff,beta),col=viridis(2)[2],lwd=2,lty=3)+
  geom_hline(yintercept = -3,lty=2)+
  geom_hline(yintercept = 3,lty=2)+
  scale_x_continuous('a-c')+
  scale_y_continuous('slope')+
  ggtitle('sd=1')


SDs <- rep(0.25,500)
D2 <- lapply(SDs,PoisMix) %>% rbindlist

D2[,beta:=6/(1+exp(-diff))-3]
g2 <- ggplot(D2,aes(diff,slope))+
  geom_point(cex=5,color=rgb(0,0,0,0.01),pch=16)+
  geom_smooth(lwd=2)+
  geom_line(aes(diff,beta),col=viridis(2)[2],lwd=2,lty=3)+
  geom_hline(yintercept = -3,lty=2)+
  geom_hline(yintercept = 3,lty=2)+
  scale_x_continuous('a-c')+
  scale_y_continuous('slope')+
  ggtitle('sd=0.25')

ggarrange(g1,g2,nrow=2)


SDs <- rep(1.5,500)
D3 <- lapply(SDs,PoisMix) %>% rbindlist

D3[,beta:=6/(1+exp(-diff))-3]
g3 <- ggplot(D3,aes(diff,slope))+
  geom_point(cex=5,color=rgb(0,0,0,0.01),pch=16)+
  geom_smooth(lwd=2)+
  geom_line(aes(diff,beta),col=viridis(2)[2],lwd=2,lty=3)+
  geom_hline(yintercept = -3,lty=2)+
  geom_hline(yintercept = 3,lty=2)+
  scale_x_continuous('a-c')+
  scale_y_continuous('slope')+
  ggtitle('sd=1.5')

save(list=ls(),file='Percolation_regression_bias_workspace')


ggarrange(g2,g1,g3,nrow=3)
ggsave('Figures/Percolation_regression_bias_Mixing.png',height=8,width=4,units='in')



set.seed(1)
n=1e2
z <- rnorm(n,sd=0.25)
xx <- seq(-10,10,by=0.1)

a=.5
b=4
c=-2.5
d=-7

Y1 <- rpois(n,lambda = rate(a,b,z))
Y2 <- rpois(n,lambda=rate(c,d,z))
Y <- Y1+Y2




png('Figures/Regression_softmax_biases_simulated2.png',width=14,height=10,units='in',res=400)
  layout(matrix(c(1,2,3,4),nrow=2,byrow=T))
  
  pcol <- rgb(0,0,0,.4)
  plot(z,Y1,pch=16,cex=2,cex.axis=2,col=pcol)
  lines(xx,rate(a,b,xx),lwd=3,lty=3)
  legend('topleft',legend = c('Observations','eta1'),pch=c(16,NA),lty=c(NA,3),col=c(pcol,'black'),lwd=3,cex=1.5)
  plot(z,Y2,pch=16,cex=2,cex.axis=2,col=pcol)
  lines(xx,rate(c,d,xx),lwd=3,lty=3)
  legend('topright',legend = c('Observations','eta2'),pch=c(16,NA),lty=c(NA,3),col=c(pcol,'black'),lwd=3,cex=1.5)
  
  gg <- glm(Y~z,family=poisson)
  ss <- segmented(gg)
  gm <- gam(Y~s(z),family=poisson)
  xx=seq(min(z)-sd(z),max(z)+sd(z),length.out = 20)
  l_true <- rate(a,b,xx)+rate(c,d,xx)
  D=data.frame('z'=xx)
  l_glm <- exp(predict(gg,newdata=D))
  l_gam <- exp(predict(gm,newdata=D))
  l_seg <- exp(predict(ss,newdata=D))
  
  cols <- viridis(4)[1:3]
  plot(z,Y,cex.axis=2,cex=2,xlim=c(min(z)-.5*sd(z),max(z)+.5*sd(z)),ylim=c(0,25),pch=16,col=rgb(0,0,0,0))
  lines(xx,l_true,lty=1)
  lines(xx,l_gam,lty=2,col=cols[2],lwd=4)
  lines(xx,l_seg,lty=3,col=cols[1],lwd=4,type='o',pch=17,cex=1.5)
  lines(xx,l_glm,lty=4,col=cols[3],lwd=4,type = 'o',pch=4,cex=1.5)
  legend(x=-0.3,y=25,legend = c('observations','true_mean','glm','gam','piecewise glm'),lwd=1.5,
         lty=c(NA,1,4,2,3),pch=c(16,NA,4,NA,17),col=c(pcol,'black',rev(cols)),cex=1.5)
  points(z,Y,pch=16,cex=2,col=pcol)
  
  
  # psi <- (a-c)/(d-b) #true switching threshold
  psi <- -0.25  ## within the dominance of Y1
  z2 <- z[z>psi]
  Y3 <- Y[z>psi]
  gg2 <- glm(Y3~z2,family=poisson)
  gm2 <- gam(Y3~s(z2),family = poisson)
  ss2 <- segmented(gg2)
  D <- data.frame('z2'=xx)
  l_gam2 <- exp(predict(gm2,newdata=D))
  l_seg2 <- exp(predict(ss2,newdata=D))
  plot(xx,l_true,type='l',lwd=2,cex=1.5,cex.axis=2,
       xlim=c(min(z)-.5*sd(z),max(z)+.5*sd(z)),ylim=c(0,25))
  points(z2,Y3,cex=2,pch=16,col=pcol)
  points(z[z<=psi],Y[z<=psi],cex=2)
  abline(v=psi)
  lines(xx,l_gam2,lty=2,col=cols[2],lwd=4)
  lines(xx,l_seg2,lty=3,col=cols[1],lwd=4,type='o',pch=17,cex=1.5)
  lines(xx,rate(coef(gg2)[1],coef(gg2)[2],xx),
        lty=4,col=cols[3],lwd=4,type = 'o',pch=4,cex=1.5)
  legend(x=-.68,y=25,legend = 'unobserved',pch=1,cex=1.5)
  legend(x=-0.2,y=25,legend = c('observations','true_mean','glm','gam','piecewise glm'),lwd=1.5,
         lty=c(NA,1,4,2,3),pch=c(16,NA,4,NA,17),col=c(pcol,'black',rev(cols)),cex=1.5)
  
dev.off()

library(gbm)






# Filtering nonlinearity --------------------------------------------------

rm(list=ls())

n=1e4
rate <- function(a,b,x) exp(a+b*x)
ilogit <- function(a,b,x) 1/(1+exp(-a-b*x))
set.seed(1)

x <- rnorm(n)
eta1 <- rate(3,-1,x)   ## this will be the Poisson canonical parameter
eta2 <- ilogit(3,0.5,x)  ## logit(p)>>0 --> Should get nonlinear, slope ~-1
eta3 <- ilogit(-3,0.5,x) ## logit(p)<<0 --> should get ~linear, slope ~-0.5


y1 <- rpois(n,lambda = eta1*eta2)
y2 <- rpois(n,lambda = eta1*eta3)


g1 <- glm(y1~x,family=poisson)
g2 <- glm(y2~x,family=poisson)
g1
g2

par(mfrow=c(1,2))
plot(x,y1,main='logit(p)>>0')
lines(sort(x),g1$fitted.values[order(x)],col='red',lwd=2)
plot(x,y2,main='logit(p)<<0')
lines(sort(x),g2$fitted.values[order(x)],col='red',lwd=2)





############# Simulate approximations
set.seed(1)
n=1e5
x <- rnorm(n,sd=0.5)

a <- 3
b <- -1
d <- 2

D <- data.table('c'=seq(-10,10,by=0.1),'slope'=0,'intercept'=0)

## at lower end, slope ~C-1. At upper end, should get slope ~-1
getCoefs <- function(cc,a=3,b=-1,d=2,x.=x){
  y <- rpois(n,lambda=rate(a,b,x)*ilogit(cc,d,x))
  gg <- glm(y~x,family=poisson)
  return(c('intercept'=coef(gg)[1],'slope'=coef(gg)[2]))
}
cl <- makeCluster(7)
clusterExport(cl,varlist=c('rate','ilogit','n','x'))
coefs <- parSapply(cl,D$c,getCoefs)
stopCluster(cl)
rm('cl')

D[,intercept:=coefs[1,]]
D[,slope:=coefs[2,]]

# for (cc in D$c){
#   y <- rpois(n,lambda=rate(a,b,x)*ilogit(cc,d,x))
#   gg <- glm(y~x,family=poisson)
#   D[c==cc]$intercept <- coef(gg)[1]
#   D[c==cc]$slope <- coef(gg)[2]
# }

cls <- viridis(3)[c(3,2,1)]
cls[2] <- 'black'

gg_intercept <- ggplot(D,aes(c,intercept))+
  geom_line(lwd=2,color=cls[2])+
  geom_abline(intercept = a,slope=1,lty=2,lwd=2,color=cls[1])+
  geom_hline(yintercept = a,lty=5,lwd=3,color=cls[3])+
  ggtitle('Intercept')+theme_gray()

gg_slope <- ggplot(D,aes(c,slope))+
  geom_line(lwd=2,color=cls[2])+
  geom_hline(yintercept = b+d,lty=2,lwd=2,color=cls[1])+
  geom_hline(yintercept = b,lty=5,lwd=3,color=cls[3])+
  ggtitle('Slope')+theme_gray()

ggarrange(gg_intercept,gg_slope,nrow = 2)
ggsave('Figures/Filtration_slope_estimation.png',height=10,width=5,units='in')


set.seed(1)
n=1e3
x <- rnorm(n)
DD <- data.table('x'=x,'y'=rpois(n,lambda= rate(a,b,x)*ilogit(2,d,x)))
g5 <- ggplot(DD,aes(x,y))+
  geom_point()+
  geom_smooth(method=glm,family=poisson,col=plasma(2)[1],show.legend=T)+
  geom_smooth(method=mgcv::gam,family=poisson,col=plasma(2)[2])+theme_bw()+
  ggtitle('beta=2')

DD <- data.table('x'=x,'y'=rpois(n,lambda= rate(a,b,x)*ilogit(-2,d,x)))
g3 <- ggplot(DD,aes(x,y))+
  geom_point()+
  geom_smooth(method=glm,family=poisson,col=plasma(2)[1],show.legend=T)+
  geom_smooth(method=mgcv::gam,family=poisson,col=plasma(2)[2])+theme_bw()+
  ggtitle('beta=-2')

DD <- data.table('x'=x,'y'=rpois(n,lambda= rate(a,b,x)*ilogit(0,d,x)))
g4 <- ggplot(DD,aes(x,y))+
  geom_point()+
  geom_smooth(method=glm,family=poisson,col=plasma(2)[1],show.legend=T)+
  geom_smooth(method=mgcv::gam,family=poisson,col=plasma(2)[2])+theme_bw()+
  ggtitle('beta=0')

library(cowplot)

ggdraw() +
  draw_plot(gg_intercept, 
            x = 0.01, y = 0.50, width = .5, height = .5) +
  draw_plot(gg_slope, 
            x = 0.01, y = .01, width = .5, height = .5) +
  draw_plot(g3, x = 0.55, y = 0.0, width = .4, height = .34) +
  draw_plot(g4, x = 0.55, y = 0.33, width = .4, height = .34) +
  draw_plot(g5, x = 0.55, y = 0.66, width = .4, height = .34) 
  
  # draw_plot(g3, x = .7, y = 0, width = .3, height = .33) +
  # draw_plot_label(label = c("Phylogeny", "A", "B","C"), size = 15,
  #                 x = c(0, 0.67, .67,.67), y = c(1, 1, 0.66,0.33))

save(list=ls(),file='Filteration_nonlinearity_demo')

dd <- expand.grid('c'=seq(-6,6,by=3),'x'=seq(-7,7,by=0.01)) %>% as.data.table
dd[,a:=a]
dd[,b:=b]
dd[,d:=d]

dd[,eta:=rate(a,b,x)*ilogit(c,d,x)]
dd[,f:=dnorm(x,sd=0.5)]

g6 <- ggplot(dd,aes(x,eta,by=rev(factor(c)),color=factor(c)))+
  geom_line(lwd=2)+
  scale_color_manual(values=viridis(length(unique(dd$c))))+
  scale_y_continuous(trans='log',breaks=10^seq(-5,5))+
  geom_abline(slope=b,intercept=a,lty=2,lwd=2)+
  theme(legend.position = c(.65,.2))


ggdraw() +
  draw_plot(gg_intercept, 
            x = 0.01, y = 0.50, width = .3, height = .5) +
  draw_plot(gg_slope, 
            x = 0.01, y = .01, width = .3, height = .5) +
  draw_plot(g6, x = 0.3, y = 0, width = .7, height = .95)

ggsave('Figures/Filtration_nonlinearities.png',height=6,width=9,units='in')
