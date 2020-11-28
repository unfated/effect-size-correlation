#packages used
if(!("devtools" %in% rownames(installed.packages()))) install.packages("devtools"); .rootdir <- "/home/unfated/Tim"; source("/home/unfated/Tim/Tmisc/inst/Startup.R")

Tmisc()
Tim.load("Rplink",lib.dir="/home/unfated/Tim")
filename <- Rfilename("sim.10000.1000", seed = "prompt")
1
library("RcppArmadillo")
library("data.table")
library("Matrix")

Tim.load("UKBBscripts", "/psychipc01/disk3/data/UKBB")
Tim.load("lassosum","/home/unfated/Tim")
Tim.load("crosspred", "/home/unfated/Tim")
#cl <- startParallel(10)

fixed.mem <- 40000

# bfile <- attachroot("~/DATA/UKBB2/data/genotype/whiteBritish1_train/whiteBritish1_train")

# bfile <- attachroot("/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train")
# ncol.bfile(bfile)
# bfile<-'/psychipc01/disk3/data/UKBB/data/genotype/imputed_0.01/maf_0.01_1'
bfile<-"/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train22"
#bfile<-'/home/unfated/chr1_block/block_1'
#bfile <- useWTCCC(.bfile = F)

N <- nrow.bfile(bfile)

set.seed(1)
ncol.bfile(bfile)
prop = c(0.01,0.125,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95)
h2=c(0.1,0.25,0.5,0.75,0.9)

# beta <- rep(0, p)
# input=c('/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train')
# impute=c('median')
# i=22
# library('KRIS')
# library('glmnet')
# library('monomvn')
# library('DMwR')
# IID=list()
# data1<-NA
data3<-NA
# bed=paste0(input,i,'.bed')
# bim=paste0(input,i,'.bim')
# fam=paste0(input,i,'.fam')

#get 22 genotype
library(data.table)
# data4 = as.data.table(data3)
# fwrite(data4,'/home/unfated/simulation/chr22_geno_20k_complete.txt')


data4=fread('/home/unfated/simulation/chr22_geno_20k_complete.txt')
# x2=read.big.matrix('/home/unfated/simulation/chr22_geno_20k_complete.txt',header=T)
data3 = as.data.frame(data4)
#x1, data2 --final genotype
#simulate phenotype
# m <- 10000
# data3 <- data2[1:m,]
x1=as.matrix(data3)

library('lsgl')
# Y <- rep(NA, N)
# Y[discovery] <- y


#phenotype done

Beta<-NA
numselect<-NA
k=10000
k2=c(1:1000)
# k2=sample(1:10000,1000)
# k2=c(1,10,100,1000)
snpid<-list()
library("lmridge") 
fulldata=data3
ncol(fulldata)
nrow(fulldata)
fulldata1 = fulldata[,which(colSums(abs(fulldata)) !=0)]
ncol(fulldata1)
useddata=fulldata1[1:k,k2]
useddata = useddata[,which(colSums(abs(useddata)) !=0)]
sum(is.na(useddata))
ncol(useddata)
nrow(useddata)
# ld=cor(useddata)-diag(1,ncol(useddata),ncol(useddata))
# length(which(ld>0.9))
useddata1=useddata
j=1
while(j<ncol(useddata1)){
  ldd=cor(useddata1[,j],useddata1[,-j])
  if(length(which(ldd>0.9))>0){
    useddata1=useddata1[,-j]}else{j=j+1}
}

ncol(useddata1)
nrow(useddata1)

sum(is.na(useddata1))
# useddata1=useddata1[,which(colSums(abs(useddata1)) !=0)]

useddata1[,709]  
useddata1=useddata1[,-709]
useddata2=useddata1[1:883,]
ncol(useddata2)
nrow(useddata2)
# var(useddata1)
#generate phenotype

for(r1 in 1:10){
  for(r2 in 1:5){
    useddata=useddata2
    p=ncol(useddata) 
    n.causal = ceiling(p*prop[r1])
    beta <- rep(0, p)
    beta[sample(p, n.causal)] <- rnorm(n.causal)
    Xb <- as.matrix(useddata)%*%beta
    herit <- h2[r2]
    error <- rnorm(length(Xb), sd=sqrt(var(Xb)/herit*(1-herit)))
    y <- Xb + error
    length(y)  
    # which(ld>0.99)
    useddata$y2<-y
    # useddata=scale(as.matrix(useddata))
    useddata=as.data.frame(useddata)
    sum(is.na(useddata[1,]))
    sum(is.na(useddata))
    useddata=na.omit(useddata)
    sum(is.na(useddata))
    # which(is.na(useddata[2,]))
    # useddata$y2
    lambdaused=(1-herit)/herit
    lambdaused
    mod <- lmridge(y2 ~ ., data = useddata, scaling = "centered", K = lambdaused) 
    Z=mod$Z
    Z=as.data.frame(Z)
    nrow(Z)
    Z=as.matrix(Z)
    nrow(Z)
    ncol(Z)
    X=as.matrix(useddata[,-ncol(useddata)])
    nrow(X)
    H=X%*%Z
    nrow(H)
    v=nrow(H)-sum(diag(H))
    # v
    ZZ=Z%*%t(Z)
    #ZZ
    head(Z%*%as.matrix(useddata[,ncol(useddata)]))
    Beta=as.vector(coef(mod))
    head(Beta)
    res=residuals(mod)
    head(res)
    #res1=scale(useddata$y2)-scale(X)%*%as.matrix(beta))
    #head(res1)
    sig=sum(res^2)/v
    se=sqrt(sig*diag(ZZ))
    t<-NA
    
    for(i in 1:length(se)){
      t[i]=coef(mod)[i+1]/se[i]
    }
    # t
    pval=2*(1-pnorm(abs(t)))
    
    write.table(pval,paste0('/home/unfated/pval0814',herit,"_",prop[r1],'.txt'))
    num=length(which(pval<0.05/length(pval)))
    num
    min(pval) 
    summary(pval)
    
    
    effect1=data.frame(beta[!is.na(pval)])
    nrow(effect1)
    pval=na.omit(pval)
    length(pval)
    # cutoff=seq(0.01)
    FPR<-NA
    FDR<-NA
    power<-NA
    FNR<-NA
    for(k in 1:10){
      BETA=as.numeric((-log10(pval)>1+0.5*k))
      # BETA
      BETA=as.data.frame(BETA)
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
      
      FPR[k]=FP/(FP+TN)
      FNR[k]=FN/(FN+TP)
      if(FP+TP==0){FDR[k]=0}
      else{FDR[k]=FP/(FP+TP)}
      power[k]=1-FN/(FN+TP)
    }
    
    perform<-data.frame(FPR,FDR,power)
    write.table(perform,paste0('/home/unfated/perform0815',herit,"_",prop[r1],'.txt'))
    # perform[[k]]=c(FPR,FNR,FDR,power,herit/n.causal,n.causal,prop[k])
    # #length(which(BETA[,1]==0))
    # #length(which(BETA[,1]!=0))
    # print(paste0('k=',k))
    # 
    # store[[l]]=perform
    # # write.table(perform,paste0('/home/unfated/simulation/powerspectrum/SGFLasso-h',l,'.txt'))
    # 
  }}



