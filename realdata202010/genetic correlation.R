library(tidyverse)
ratio.result<-list()
for(k in 1:5000){
a.all<-seq(0,3,by=0.25)
b.all<-seq(5,30,by=1)

s=c(sample(seq(-3,3,by=0.5),1),sample(seq(-3,3,by=0.5),1))#s1,s2>=0
ss=c(mean(s),mean(s),s[1]-1,s[2]-1,mean(s)-1,mean(s)-1,s[1],s[2])
aa=sample(a.all,1)
bb=sample(b.all,1)
while (aa+s[1]<=1|aa+s[2]<=1){aa=aa+0.25}
para=c(aa,bb) #a<=2, b>=5
nom=1
denom=1
for(i in 1:4){
  a=para+ss[i]
  nom=nom*beta(a[1],a[2])
  b=para+ss[i+4]
  denom=denom*beta(b[1],b[2])
}
ratio.result[[k]]=c(para,s,aa/(aa+bb),nom/denom)
}
ratio.result1<-ratio.result %>% as.data.frame(.) %>% t() %>% data.frame()
colnames(ratio.result1)=c('alpha','beta','S1','S2','mean.MAF','ratio')
rownames(ratio.result1)=c(1:nrow(ratio.result1))

par(mfrow=c(1,1))
plot(ratio.result1$S1,ratio.result1$ratio)
plot(ratio.result1$S1,ratio.result1$ratio,xlab='s1 (s2)',ylab='ratio')
persp(ratio.result1$S1,ratio.result1$S2,ratio.result1$ratio, col='blue')
install.packages('plot3D')
library(plot3D)
install.packages('scatterplot3d')
library(scatterplot3d)
scatterplot3d(ratio.result1$S1,ratio.result1$S2,ratio.result1$ratio,highlight.3d = T,main='ratio vs (s1,s2)',xlab='s1',ylab='s2',zlab='ratio')
library(car)
symbols(ratio.result1$S1,ratio.result1$S2,circle=ratio.result1$ratio,main='ratio vs (s1,s2)',xlab='s1',ylab='s2',inches=0.2)
#simulation----
h1=0.8
rho_b=0.7
beta<-mvrnorm(1000,c(0,0),matrix(c(h1,h1*rho_b,h1*rho_b,h1),ncol=2))
MAF=rbeta(1000,1.1,10)
beta2<-beta/sqrt(MAF*(1-MAF))
cov(beta2)
cor(beta2)
