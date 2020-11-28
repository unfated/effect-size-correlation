# impute as 2
ccc='/home/unfated/Tim/'
itr=c(1:10)
option='lasso'
Lasso(itr,ccc)

# BLasso

# impution compare
list<-paste0('/home/unfated/Tim/beta.ls_fit',c(1:133),'.txt')
BETA1=c()
for(i in 1:10){
  beta<-read.delim(list[i], stringsAsFactors=FALSE,sep = ' ')
  beta=beta[-1,]
  BETA1=c(BETA1,beta)}

list<-paste0('/home/unfated/Tim/beta2.ls_fit',c(1:133),'.txt')
BETA2=c()
for(i in 1:10){
  beta<-read.delim(list[i], stringsAsFactors=FALSE,sep = ' ')
  beta=beta[-1,]
  BETA2=c(BETA2,beta)}

length(BETA1)
length(BETA2)
which(BETA1!=BETA2)
which(BETA2!=BETA1)
write.table(trueB,'/home/unfated/Tim/effectsize_chr1_remove_naSNP.txt')

trueb=read.table('/home/unfated/Tim/effectsize_chr1_remove_naSNP.txt')
head(trueb)
trueB<-trueb[,1]
holder1='/home/unfated/Tim/beta.ls_fit'
holder2='/home/unfated/Tim/beta2.ls_fit'
Rate(trueB=trueB,holder=holder1,num=10)
Rate(trueB=trueB,holder=holder2,num=10)
head(BETA1)
head(BETA2)
