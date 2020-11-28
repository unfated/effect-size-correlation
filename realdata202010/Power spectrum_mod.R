# power spectrum
#packages used
if(!("devtools" %in% rownames(installed.packages()))) install.packages("devtools"); .rootdir <- "/home/unfated/Tim"; source("/home/unfated/Tim/Tmisc/inst/Startup.R")
#17
## fixed.mem = 40gb
# a package
Tmisc()
Tim.load("Rplink",lib.dir="/home/unfated/Tim")
filename <- Rfilename("sim.10000.1000", seed = "prompt")
1
#0
# test.n <- 2000
# validation.n <- 2000
# m <- as.integer(nameoption(filename, 2))
# n.causal <- as.integer(nameoption(filename, 3))
library("RcppArmadillo")
library("data.table")
library("Matrix")

Tim.load("UKBBscripts", "/psychipc01/disk3/data/UKBB")
Tim.load("lassosum","/home/unfated/Tim")
Tim.load("crosspred", "/home/unfated/Tim")
#cl <- startParallel(10)

fixed.mem <- 40000
#bfile <- attachroot("~/DATA/UKBB2/data/genotype/whiteBritish1_train/whiteBritish1_train")
#bfile <- attachroot("/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train")
# bfile<-'/psychipc01/disk3/data/UKBB/data/genotype/imputed_0.01/maf_0.01_1'
#bfile<-"/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train1"
bfile<-'/home/unfated/chr1_block/block_1'
#bfile <- useWTCCC(.bfile = F)

N <- nrow.bfile(bfile)
m <- 1000
set.seed(1)
discovery.sample <- sample(N, m)
# others.sample <- setdiff(1:N, discovery.sample)
# validation.sample <- sample(others.sample, validation.n)
# others.sample <- setdiff(1:N, c(discovery.sample, validation.sample))
# 
discovery <- logical.vector(discovery.sample, N)
# validation <- logical.vector(validation.sample, N)

p <- ncol.bfile(bfile)
prop = c(0.01,0.125,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95)

# beta <- rep(0, p)

option=c('lasso')
holder=c('/home/unfated/Tim/')
holder2=c('/home/unfated/simulation/')
input=c('/home/unfated/chr1_block/block_')
impute=c('median')
i=1
library('KRIS')
library('glmnet')
library('monomvn')
library('DMwR')
IID=list()
data1<-NA
data3<-NA
bed=paste0(input,i,'.bed')
bim=paste0(input,i,'.bim')
fam=paste0(input,i,'.fam')
block=read.bed(bed,bim,fam)
z_1=block$snp[discovery,]
data1=data.frame(z_1)
id=c()
for(k in 1:ncol(data1)){
  if(sum(is.na(data1[,k]))>nrow(data1)/2)
    id=c(id,k)
}
IID[[i]]=id
data3=data1
if(length(id>0))
  data3=data1[,-id]
if(impute=='1'){
  for(j in 1:ncol(data3)){
    data3[is.na(data3[,j]),j]=1
  }}
if(impute=='2'){
  for(j in 1:ncol(data3)){
    data3[is.na(data3[,j]),j]=2
  }}
if(impute=='mean'){
  for(j in 1:ncol(data3)){
    me=mean(data3[,j],na.rm=T)
    data3[is.na(data3[,j]),j]=me
  }}
if(impute=='median'){
  for(j in 1:ncol(data3)){
    md=median(data3[,j],na.rm=T)
    data3[is.na(data3[,j]),j]=md
  }}
if(impute=='0'){
  for(j in 1:ncol(data3)){
    data3[is.na(data3[,j]),j]=0
  }}
if(impute=='KNN_avg'){
  data2=knnImputation(data3, k = 5, scale = T, meth = "weighAvg", distData = NULL)
  data3=data2}
if(impute=='KNN_med'){
  data2=knnImputation(data3, k = 5, scale = T, meth = "median", distData = NULL)
  data3=data2}
x1=as.matrix(data3)




perform=list() 
p=ncol(data3) 
h1=0.001
h2=0.25
for(k in 1:10){
  n.causal = ceiling(p*prop[k])
  
  library('STAT')
  mu=rep(0,n.causal)
  # set.seed(1)
  # V=rWishart(1,n.causal,diag(n.causal))
  # V2=rWishart(1,n.causal,V[,,1])
  #corr=cov2cor(as.matrix(V2[,,1]))
  corr=diag(n.causal)
  library(MASS)
  #cov=matrix(c(1,.5,.15,.15,0,.5,1,0.15,0.15,0,.15,.15,1,.25,0,.15,.15,.25,1,0,0,0,0,0,1),nrow=5,ncol=5)
  set.seed(1)
  traitnum=100
  effect=mvrnorm(n = traitnum, mu=mu, Sigma=corr, tol = 1e-6, empirical = FALSE)
  # write.table(effect,'/home/unfated/simulation/effect.txt')
  # colnames(corr)<-paste0('SNP',1:ncol(corr))
  # rownames(corr)<-paste0('SNP',1:ncol(corr))
  # write.table(corr,'/home/unfated/simulation/correlation.txt')
  
  
  

  Betalasso=list()
  d=n.causal
  ccc=sample(c(1:p),d)
  for(j in 1:traitnum){
    beta <- rep(0, p)
    
    beta[ccc]<- effect[j,]
    Xb<-x1%*%beta
    herit <- runif(1,min=h1,max=h2)
    error <- rnorm(m, sd=sqrt(var(Xb)/herit*(1-herit)))
    y <- Xb + error
    cv1 = cv.glmnet(x1, y)
    beta.ls_fit1=coef(cv1, s = "lambda.1se")
    write.table(as.matrix(beta.ls_fit1),paste0(holder2,'beta_',impute,'.ls_fit',i,'_trait',j,'.txt'),col.names = F)
    Betalasso[[j]]=data.frame(as.matrix(beta.ls_fit1))[,1][-1]
    print(j)}
  
  #write.table(BETA,'/home/unfated/simulation/lassobeta.txt')
  
  l=nrow(effect)*ncol(effect)
  BETA=t(as.data.frame(Betalasso))
  effect1=data.frame(matrix(0,nrow=traitnum,ncol=p))
  effect1[ccc]=effect
  TN=0
  TP=0
  FP=0
  FN=0
  for(i in 1:nrow(effect1)){
    for(j in 1:ncol(effect1)){
      if(effect1[i,j]==0&& BETA[i,j]==0)
        TN=TN+1
      if(effect1[i,j]==0&& BETA[i,j]!=0)
        FP=FP+1
      if(effect1[i,j]!=0&& BETA[i,j]==0)
        FN=FN+1
      if(effect1[i,j]!=0&& BETA[i,j]!=0)
        TP=TP+1
    }}
  print(paste0('FPR=',FP/(FP+TN)))
  print(paste0('FNR=',FN/(FN+TP)))
  print(paste0('FDR=',FP/(FP+TP)))
  print(paste0('Power=',1-FN/(FN+TP)))
  
  FPR=FP/(FP+TN)
  FNR=FN/(FN+TP)
  FDR=FP/(FP+TP)
  power=1-FNR
  perform[[k]]=c(FPR,FNR,FDR,power,(h1+h2)/(2*n.causal),n.causal,prop[k])
  #length(which(BETA[,1]==0))
  #length(which(BETA[,1]!=0))
  print(paste0('k=',k))
}
write.table(perform,'/home/unfated/simulation/powerspectrum/uniLasso-h1.txt')