library(tidyverse)

prop_sort <- function(type) {
  # type='sim'
  data=list()
  for(i in 1:10){
  data[[i]] <- read.csv(paste0('C:/Users/22367/Desktop/thesis/fine mapping real data/20201112 result/prop.logcor20201113.',type,i,".csv"),header= T,row.names = 1)
  colnames(data[[i]])=c(10,20,50,100,200,300,400,500)}
  
  prop.mean=data[[1]]
  prop.sd=data[[1]]
  
  for(i in 1:nrow(prop.mean)){
    for(j in 1:ncol(prop.mean)){
      a=c()
  for(k in 1:10){
    a[k]=data[[k]][i,j]
  }
  prop.mean[i,j]=mean(a)
  prop.sd[i,j]=sd(a)}}
  

  b=rownames(data[[1]])
  
  aa=gather(prop.mean,key='sample_size',value='frequency')
  bb=gather(prop.sd,key='sample_size',value='sd')
  bb
  aa$bin=rep(b,8)
  aa$sd=bb$sd
  aa$bin=factor(aa$bin)
  aa$sample_size=factor(aa$sample_size)
return(aa)

}

final.data.real=prop_sort('real')
final.data.sim=prop_sort('sim')
final.data.perm=prop_sort('perm')

final.data.perm$group=paste0('permu: ',final.data.perm$sample_size)
final.data.real$group=paste0('real: ',final.data.perm$sample_size)
final.data.sim$group=paste0('simu: ',final.data.perm$sample_size)
rbind(final.data.perm,final.data.real,final.data.sim) %>% 
ggplot(., aes(x=bin,y=frequency,fill=group)) + theme_bw()+
  geom_col(color='black',position='dodge')
final.data=rbind(subset(final.data.perm,group=='permu: 500'),final.data.real,subset(final.data.sim,group=='simu: 500')) 

ggplot(final.data, aes(x=bin,y=frequency,fill=group)) + theme_bw()+
  geom_col(color='black',position='dodge')

a=sort(c(10, 20,50,100,200,300,400,500),decreasing = T)
a
b=c('permu: 500','simu: 500',paste0('real: ',a))
b
c=unique(final.data$bin)
final.data$group=factor(final.data$group,levels=b)
c
final.data$bin=factor(final.data$bin,levels=c)
pd<- position_dodge(0.9)
ggplot(subset(final.data,bin!='(-1,-0.5]'&bin!='(-0.5,-0.25]'), aes(x=bin,y=frequency,fill=group)) + theme_bw()+
  geom_col(color='black',position='dodge')+
  geom_errorbar(aes(ymin = frequency-sd, ymax = frequency + sd),position=pd,width=0.2)+
  scale_fill_discrete(name = 'Distribution groups',labels=c('Null by permutation: n=500','Null by simulation: Pval~uniform(0,1)',paste0('real data: n=',a)))+
     theme(                
    legend.justification=c(1,0),
        legend.position=c(0.2,0.5))  +
  xlab('bins: Pearson sample correlation r')+
  ylab(expression('frequency:'~'n'^2~'correlations of n'^2~'gene pairs'))
# a=runif(1000,-0,1) %>% 
# -log(.,10)
# b=runif(1000,-0,1) %>% 
#   -log(.,10)
# cor(a,b)
# cor(runif(1000,-0,1),runif(1000,-0,1))
