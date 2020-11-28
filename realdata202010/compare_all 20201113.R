#packages used----
if(!("devtools" %in% rownames(installed.packages()))) install.packages("devtools"); .rootdir <- "/home/unfated/Tim"; source("/home/unfated/Tim/Tmisc/inst/Startup.R")

Tmisc()
Tim.load("Rplink",lib.dir="/home/unfated/Tim")
filename <- Rfilename("sim.10000.1000", seed = "prompt")
1
library("RcppArmadillo")
library("data.table")
library("Matrix")
library('tidyverse')
Tim.load("UKBBscripts", "/psychipc01/disk3/data/UKBB")
Tim.load("lassosum","/home/unfated/Tim")
Tim.load("crosspred", "/home/unfated/Tim")

library('lsgl')
library("lmridge")
library(data.table)
library(devtools)
library(glmnet)
#devtools::install_github("nielsrhansen/sglOptim", build_vignettes = TRUE)
#devtools::install_github("nielsrhansen/lsgl", build_vignettes = TRUE)
#load genotype----

bfile<-"/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train22"

nrow.bfile(bfile)

set.seed(1)
ncol.bfile(bfile)
#get chr 22 genotype


data4=fread('/home/unfated/simulation/chr22_geno_20k_complete.txt')
data3 = as.data.frame(data4)


#individuals number
k1=c(1:10000)
k2=c(1:1000)
 
extract_geno<-function(data3,k1,k2,n.ind,n.snp){
fulldata=data3
ncol(fulldata)
nrow(fulldata)
fulldata1 = fulldata[,which(colSums(abs(fulldata)) !=0)]
ncol(fulldata1)
useddata=fulldata1[k1,k2]
#first
useddata = useddata[,which(colSums(abs(useddata)) !=0)]
sum(is.na(useddata))
ncol(useddata)
nrow(useddata)
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
var.data=apply(useddata1,2,var) 

useddata1=useddata1[,-which(var.data==0)]
useddata2=useddata1[n.ind,n.snp]
return(useddata2)}


# useddata2: data for ridge





#generate phenotype----


generate_pheno<-function(useddata2,r1,r2){
  Beta<-NA
  numselect<-NA
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
return(list(y,beta))}

#method: ridge----
ridge.regress<-function(useddata,y,lambdaused){
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
  
    
    # lambdaused
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
    
return(list(t,pval))}


#eval function
perform.eval<-function(beta,pval,cutoff){ 
    effect1=data.frame(beta[!is.na(pval)])
    nrow(effect1)
    pval=na.omit(pval)
    length(pval)
    # cutoff=seq(0.01)
    FPR<-NA
    FDR<-NA
    power<-NA
    FNR<-NA

      BETA=as.numeric((-log10(pval)>cutoff))
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
     # print(paste0('FPR=',FP/(FP+TN)))
      #print(paste0('FNR=',FN/(FN+TP)))
      #print(paste0('FDR=',FP/(FP+TP)))
      #print(paste0('Power=',1-FN/(FN+TP)))
      
      FPR=FP/(FP+TN)
      FNR=FN/(FN+TP)
      if(FP+TP==0){FDR=0}
      else{FDR=FP/(FP+TP)}
      power=1-FN/(FN+TP)
      F1.score=2*power*(1-FDR)/(power+(1-FDR))

    
    perform<-c(sum(BETA),cutoff,FPR,FDR,power,F1.score)
    return(perform)}
#method gwas----

#gwas

lm_getp<-function(x){a=summary(lm(y~x))$coefficients
b=1
if(nrow(a)>=2){b=a[2,4]}
return(b)}




cutoff.perm<-function(ridge.pval) {
perform.ridge=matrix(0,nrow=100,ncol=6)
for(i in 1:100){
a=seq(0.1,10,length.out=100)[i]
perform.ridge[i,]=perform.eval(beta,ridge.pval,a)
}
perform.ridge=data.frame(perform.ridge)
colnames(perform.ridge)=c('No.discovery','cutoff','FPR','FDR','power','F1.score')
best=subset(perform.ridge,perform.ridge[,6]==max(perform.ridge[,6]))
return(list(perform.ridge,best))}


# conduct----
k1=c(1:10000)
k2=c(1:1400)
n.ind=c(1:1000)
n.snp=c(1:1000)
useddata2=extract_geno(data3,k1,k2,n.ind,n.snp)
ncol(useddata2)

prop = c(0.01,0.0125,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95)
h2=c(0.1,0.25,0.5,0.75,0.9)





#univariate lasso method----
cutoff.perm2<-function(lasso.fit) {

perform.ridge=matrix(0,nrow=length(lasso.fit$lambda),ncol=6)
for(i in 1:length(lasso.fit$lambda)){
beta.lasso=as.matrix(coef(lasso.fit, s = lasso.fit$lambda[i]))[,1][-1]
perform.ridge[i,]=c(lasso.fit$lambda[i],perform.eval2(beta,beta.lasso))
}
perform.ridge=data.frame(perform.ridge)
colnames(perform.ridge)=c('lambda','No.discovery','FPR','FDR','power','F1.score')
best=subset(perform.ridge,perform.ridge[,6]==max(perform.ridge[,6]))
return(list(perform.ridge,best))}

perform.eval2<-function(beta,pval){ 

    effect1=data.frame(beta[!is.na(pval)])
    #nrow(effect1)
    pval=na.omit(pval)
    #length(pval)
    # cutoff=seq(0.01)
    FPR<-NA
    FDR<-NA
    power<-NA
    FNR<-NA

      BETA=as.numeric(pval!=0)
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
     # print(paste0('FPR=',FP/(FP+TN)))
      #print(paste0('FNR=',FN/(FN+TP)))
      #print(paste0('FDR=',FP/(FP+TP)))
      #print(paste0('Power=',1-FN/(FN+TP)))
      
      FPR=FP/(FP+TN)
      FNR=FN/(FN+TP)
      if(FP+TP==0){FDR=0}
      else{FDR=FP/(FP+TP)}
      power=1-FN/(FN+TP)
      F1.score=2*power*(1-FDR)/(power+(1-FDR))

    
    perform<-c(sum(BETA),FPR,FDR,power,F1.score)
    return(perform)}

#10 replication lasso vs ridge vs gwas

result.all1<-list()
result.all2<-list()
result.all3<-list()

(1,2)
(1,4)
(5,2)
(5,4)


para=list(c(1,2),c(1,4),c(5,2),c(5,4))

for(rr in 1:length(para)){

#prop
r1=para[[rr]][1]
#herit
r2=para[[rr]][2]
print(c('rr=',rr))
print(c('pi=',prop[r1]))
print(c('h2=',h2[r2]))
perform.lasso.rep=data.frame(matrix(0,nrow=10,ncol=6))
colnames(perform.ridge.rep)=c('lambda','No.discovery','FPR','FDR','power','F1.score')

perform.ridge.rep=data.frame(matrix(0,nrow=10,ncol=6))
colnames(perform.ridge.rep)=c('No.discovery','cutoff','FPR','FDR','power','F1.score')
perform.gwas.rep=data.frame(matrix(0,nrow=10,ncol=6))
colnames(perform.gwas.rep)=c('No.discovery','cutoff','FPR','FDR','power','F1.score')

for(ss in 1:10){
print(c('ss=',ss))
phenotype=generate_pheno(useddata2,r1,r2)
y=phenotype[[1]]
beta=phenotype[[2]]

lambda=(1-h2[r2])/h2[r2]*length(n.snp)
ridge.result=ridge.regress(useddata2,y,lambda)
ridge.pval=ridge.result[[2]]
ridge.pval[is.na(ridge.pval)]=1


gwas.pval=apply(useddata2,2,lm_getp)

perform.gwas<-cutoff.perm(gwas.pval)
perform.ridge<-cutoff.perm(ridge.pval)
#ridge and gwas
perform.gwas.rep[ss,]=perform.gwas[[2]][1,]
perform.ridge.rep[ss,]=perform.ridge[[2]][1,]

#lasso
nlambda=50
  lasso.fit= glmnet(as.matrix(useddata2), y,alpha=1,nlambda=nlambda)

 perform.lasso<-cutoff.perm2(lasso.fit)
perform.lasso.rep[ss,]= perform.lasso[[2]][1,]}

print(c('lasso',apply(perform.lasso.rep,2,mean)))
print(c('ridge',apply(perform.ridge.rep,2,mean)))
print(c('gwas',apply(perform.gwas.rep,2,mean)))
result.all1[[rr]]=perform.ridge.rep
result.all2[[rr]]=perform.gwas.rep
result.all3[[rr]]=perform.lasso.rep
}

write.table(result.all1,'/home/unfated/realdata/simulation20201113/ridge.10rep.20201113.txt')
write.table(result.all2,'/home/unfated/realdata/simulation20201113/gwas.10rep.20201113.txt')
write.table(result.all3,'/home/unfated/realdata/simulation20201113/lasso.10rep.20201113.txt')

perform.result.all=data.frame(
  apply(result.all2[[1]],2,mean)[-2],
  apply(result.all1[[1]],2,mean)[-2],
  apply(result.all3[[1]],2,mean)[-1],
 apply(result.all2[[3]],2,mean)[-2],
  apply(result.all1[[3]],2,mean)[-2],
  apply(result.all3[[3]],2,mean)[-1],
  apply(result.all2[[2]],2,mean)[-2],
  apply(result.all1[[2]],2,mean)[-2],
  apply(result.all3[[2]],2,mean)[-1],
  apply(result.all2[[4]],2,mean)[-2],
  apply(result.all1[[4]],2,mean)[-2],
  apply(result.all3[[4]],2,mean)[-1])
  colnames(perform.result.all)=rep(c('LR+T','Ridge+T','Lasso'),4)
  write.csv(round(perform.result.all,2),'/home/unfated/realdata/simulation20201113/perform.compare.all3.20201113.csv')