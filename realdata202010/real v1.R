# real data
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
#bfile <- attachroot("~/DATA/UKBB2/data/genotype/whiteBritish1_train/whiteBritish1_train")
#bfile <- attachroot("/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train")
# bfile<-'/psychipc01/disk3/data/UKBB/data/genotype/imputed_0.01/maf_0.01_1'
bfile<-"/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train22"
#bfile<-'/home/unfated/chr1_block/block_1'
#bfile <- useWTCCC(.bfile = F)

N <- nrow.bfile(bfile)
m <- 200000
set.seed(1)
discovery.sample <- sample(N, m)
discovery <- logical.vector(discovery.sample, N)
ncol.bfile(bfile)
prop = c(0.01,0.125,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95)

# beta <- rep(0, p)
input=c('/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train')
impute=c('median')
i=22
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
# block=read.bed(bed,bim,fam)
# z_1=block$snp[discovery,]
# data1=data.frame(z_1)
# id=c()
# for(k in 1:ncol(data1)){
#   if(sum(is.na(data1[,k]))>nrow(data1)/2)
#     id=c(id,k)
# }
# IID[[i]]=id
# data3=data1
# if(length(id>0))
#   data3=data1[,-id]
# if(impute=='1'){
#   for(j in 1:ncol(data3)){
#     data3[is.na(data3[,j]),j]=1
#   }}
# if(impute=='2'){
#   for(j in 1:ncol(data3)){
#     data3[is.na(data3[,j]),j]=2
#   }}
# if(impute=='mean'){
#   for(j in 1:ncol(data3)){
#     me=mean(data3[,j],na.rm=T)
#     data3[is.na(data3[,j]),j]=me
#   }}
# if(impute=='median'){
#   for(j in 1:ncol(data3)){
#     md=median(data3[,j],na.rm=T)
#     data3[is.na(data3[,j]),j]=md
#   }}
# if(impute=='0'){
#   for(j in 1:ncol(data3)){
#     data3[is.na(data3[,j]),j]=0
#   }}
# if(impute=='KNN_avg'){
#   data2=knnImputation(data3, k = 5, scale = T, meth = "weighAvg", distData = NULL)
#   data3=data2}
# if(impute=='KNN_med'){
#   data2=knnImputation(data3, k = 5, scale = T, meth = "median", distData = NULL)
#   data3=data2}
# 
# 
# sum(is.na(data3))
#get 22 genotype
library(data.table)
# data4 = as.data.table(data3)
# fwrite(data4,'/home/unfated/simulation/chr22_geno_20k_complete.txt')


data4=fread('/home/unfated/simulation/chr22_geno_20k_complete.txt')
# x2=read.big.matrix('/home/unfated/simulation/chr22_geno_20k_complete.txt',header=T)
data3 = as.data.frame(data4)
x1=as.matrix(data3)
library('lsgl')
#trait
traitnum = 10
p=ncol(data3) 
YY=matrix(0,nrow=m,ncol=traitnum)
pheno=read.delim("/psychipc01/disk3/data/UKBB/data/phenotype/UKBBphenos.txt")
IDlist=read.table(fam)
head(pheno$FID)
head(IDlist$V1)
colnames(pheno)
sum(which(pheno$FID %in% IDlist$V1))
newpheno=pheno[which(pheno$FID %in% IDlist$V1),]
colnames(newpheno)
head(newpheno)
nrow(newpheno)
YY1=newpheno[discovery,]
nrow(YY1)
nrow(YY2)

YY=YY1[,c(3,4,8,9,15:20)]
YY2=newpheno[,c(3,4,8,9,15:20)]
colnames(YY2)
colnames(YY)
    for(j in 1:traitnum){
      md=median(YY[,j],na.rm=T)
      YY[is.na(YY[,j]),j]=md}

# YY for lasso

for(j in 1:traitnum){
  md=median(YY2[,j],na.rm=T)
  YY2[is.na(YY2[,j]),j]=md}
for(j in 5:10){
  YY2[,j]=2+YY2[,j]}

sum(is.na(YY))
sum(is.na(YY2))

#YY2 for GWAS
# write.table(data3,'/home/unfated/simulation/chr22_geno_20k.txt')
#clumping

clump <- function(pvals, genotype, pos, SNP=paste0("SNP", 1:length(pvals)), 
                  r2=0.2, window.kb=250000, p1=1, p2=1, trace=0) {
  #' Function to mimic Plink's --clump command
  #' It's not equivalent to plink's command. 
  
  P <- length(pvals)
  stopifnot(P == ncol(genotype))
  stopifnot(P == length(pos))
  stopifnot(P == length(SNP))
  stopifnot(all(sort(pos) == pos))
  stopifnot(!any(is.na(pos)))
  sd0 <- apply(genotype,MARGIN = 2,function(x) sd(x,na.rm = T)) == 0
  tickedoff <- is.na(pvals) | sd0
  pvals[tickedoff] <- Inf
  
  order <- order(pvals)
  results <- data.frame(SNP=SNP, p=NA, clumped.SNPs = "")
  
  for(i in 1:P) {
    snp=order[i]
    if(tickedoff[snp]) next
    pval <- pvals[snp]
    # if(trace > 0) {
    #   message <- cat(paste0("i=", i, ", pval=", pval, ", %marked=",
    #                         sum(tickedoff)/length(tickedoff), "\n"))
    # }
    if(pval > p1) break
    po <- pos[snp]
    toinclude <- pos >= po - window.kb * 1000 & pos <= po + window.kb * 1000 &
      pvals < p2 & pos != pos[snp] & !tickedoff
    R2 <- cor(genotype[,snp], genotype[,toinclude], 
              use = "pairwise.complete.obs")^2
    R2[is.na(R2)] <- 0
    toinclude[toinclude] <- R2 >= r2
    
    results$SNP[i] <- SNP[snp]
    results$clumped.SNPs[i] <- paste(SNP[toinclude], collapse = ",")
    results$p[i] <- pvals[snp]
    
    tickedoff[toinclude] <- TRUE
    tickedoff[snp] <- TRUE
  }
  return(results[!is.na(results$p), ])
}

#GWAS
#perSNP GWAS
numGWAS<-NA
numClump<-NA
SNPinfo<-list()
clumpSNPinfo<-list()

for(j in 1:traitnum){
  
  Y=YY2[,j]
  sum <- plink(bfile=bfile, cmd="--glm", plink2=T, 
               pheno=Y, silent=F, keep=discovery
  )
  summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)
  SNP<-summ[,c('ID','P','POS')]
  SNP$P[is.na(SNP$P)] <- 1
  newSNP = subset(SNP,P < 10^(-6))
  # SNP_POS = as.numeric(rownames(newSNP))
  # genotype=data1[,SNP_POS]
  pvals=newSNP$P
  pos=newSNP$POS
   numGWAS[j]=sum(pvals<10^(-7))
# if(sum(pvals<5*10^(-8))>0){
#   clump_p=clump(pvals,genotype,pos,r2=0.5)
#   pval_new=clump_p$p
#   numClump[j]=sum(pval_new<5*10^(-8))
#   clumpSNPinfo[[j]]=pval_new[which(pval_new<5*10^(-8))]}
  SNPinfo[[j]]=subset(SNP,P < 5*10^(-8))
  
}
numGWAS
numClump
write.table(as.character(SNPinfo),'/home/unfated/simulation/SNPinfo.txt')
#pure Multivariate Lasso
y1=as.matrix(YY[,1])
    lambda<-lsgl::lambda(x1,y1, alpha=0, lambda.min=0.05)
    fit<-lsgl::fit(x = x1[1:10000,], y = y1, alpha = 1, lambda = 0.1)
    fit2<-lsgl::cv(x = x1[1:10000,], y = y1, alpha = 1, lambda = 0.1)
    betalasso=as.matrix(coef(fit, traitnum)[,-1])
    
#univariate lasso
    #gaussian
    Beta<-list()
    numselect<-NA
    k=20000
    snpid<-list()
    for(i in c(1:6,10)){
      y2=YY[,i]
      fit_cv = cv.glmnet(x1[1:k,],as.matrix(y2)[1:k,],family = 'gaussian',dfmax = ceiling(ncol(x1)*0.1),nfolds = 5,nlambda = 200)
      lambda0=fit_cv$lambda.1se
      beta_cv=coef(fit_cv,s='lambda.1se')
      print(sum(beta_cv!=0))
      print(lambda0)
      fit_lar = glmnet(x1[1:10*k,],as.matrix(y2)[1:10*k,],lambda=lambda0,family = 'gaussian',dfmax = ceiling(ncol(x1)*0.1))
      beta = coef(fit_lar)
      print(sum(beta!=0))
      print(which(beta!=0))
      snpid[[i]]=which(beta!=0)
      Beta[[i]]=as.data.frame(as.matrix(beta))
      write.table(Beta[[i]],paste0('/home/unfated/simulation/coef_lar_trait',i,'.txt'))
    }
  write.table(as.character(snpid),'/home/unfated/simulation/snpid_lar.txt')
  write.table(as.character(Beta),'/home/unfated/simulation/Beta_lar.txt')
    
    #binomial poisson
  # for(i in 10:10){
  #   y2=YY[,i]
  #   fit_cv = cv.glmnet(x1[1:k,],as.factor(y2)[1:k],family = 'binomial',dfmax = ceiling(ncol(x1)*0.1),nfolds = 5,nlambda = 200)
  #  lambda0=fit_cv$lambda.1se
  #    beta_cv=coef(fit_cv,s='lambda.1se')
  #   print(sum(beta_cv!=0))
  #   print(lambda0)
  #   fit_lar = glmnet(x1[1:10*k,],as.factor(y2)[1:10*k],lambda=lambda0,family = 'binomial',dfmax = ceiling(ncol(x1)*0.1))
  #   beta = coef(fit_lar)
  #   print(sum(beta!=0))
  #   print(which(beta!=0))
  #   snpid[[i]]=which(beta!=0)
  #   Beta[[i]]=as.data.frame(as.matrix(beta))
  #   write.table(Beta[[i]],paste0('/home/unfated/simulation/coef_lar_trait',i,'.txt'))
  # }
  #   YY4=YY+0.01
  #   for(i in 5:10){
  #     y2=YY4[,i]
  #     fit_cv = cv.glmnet(x1[1:k,],as.matrix(y2)[1:k,],family = 'gaussian',dfmax = ceiling(ncol(x1)*0.1),nfolds = 5,nlambda = 200)
  #     lambda0=fit_cv$lambda.1se
  #     beta_cv=coef(fit_cv,s='lambda.1se')
  #     print(sum(beta_cv!=0))
  #     print(lambda0)
  #     fit_lar = glmnet(x1[1:10*k,],as.matrix(y2)[1:10*k,],lambda=lambda0,family = 'gaussian',dfmax = ceiling(ncol(x1)*0.1))
  #     beta = coef(fit_lar)
  #     print(sum(beta!=0))
  #     print(which(beta!=0))
  #     snpid[[i]]=which(beta!=0)
  #     Beta[[i]]=as.data.frame(as.matrix(beta))
  #     write.table(Beta[[i]],paste0('/home/unfated/simulation/coef_lar_trait',i,'.txt'))
  #   }
   
  #el net
  #gaussian
  Beta<-list()
  numselect<-NA
  k=20000
  snpid_lar=snpid
  snpid<-list()
  for(i in c(6,10)){
    y2=YY[,i]
    fit_cv = cv.glmnet(x1[1:k,],as.matrix(y2)[1:k,],family = 'gaussian',dfmax = ceiling(ncol(x1)*0.1),nfolds = 5,nlambda = 200,alpha = 0.5)
    lambda0=fit_cv$lambda.1se
    beta_cv=coef(fit_cv,s='lambda.1se')
    print(sum(beta_cv!=0))
    print(lambda0)
    fit_lar = glmnet(x1[1:10*k,],as.matrix(y2)[1:10*k,],lambda=lambda0,family = 'gaussian',dfmax = ceiling(ncol(x1)*0.1),alpha=0.5)
    beta = coef(fit_lar)
    print(sum(beta!=0))
    print(which(beta!=0))
    snpid[[i]]=which(beta!=0)
    Beta[[i]]=as.data.frame(as.matrix(beta))
    write.table(Beta[[i]],paste0('/home/unfated/simulation/coef_eln_trait',i,'.txt'))
  }
  snpid_eln=snpid
  write.table(as.character(snpid),'/home/unfated/simulation/snpid_eln.txt')
  write.table(as.character(Beta),'/home/unfated/simulation/Beta_eln.txt')
  
  
  #ridge + thresholding
  #gaussian
  Beta<-list()
  numselect<-NA
  k=20000
  snpid<-list()
  for(i in c(1:6,10)){
  
    y2=YY[,i]
    # fit_cv = cv.glmnet(x1[1:k,],as.matrix(y2)[1:k,],family = 'gaussian',nfolds = 5,nlambda = 200,alpha = 0)
    # lambda0=fit_cv$lambda.1se
    # beta_cv=coef(fit_cv,s='lambda.1se')
    # print(sum(beta_cv!=0))
    # print(lambda0)
    lambda0=0.2
    fit_lar = glmnet(x1[1:10*k,],as.matrix(y2)[1:10*k,],lambda=lambda0,family = 'gaussian',alpha=0)
    beta = coef(fit_lar)
    print(sum(beta!=0))
    print(sum(abs(beta)>0.01)-1)
    print(sum(abs(beta)>0.003333)-1)
    print(sum(abs(beta)>0.001)-1)
    #print(which(beta!=0))
    snpid[[i]]=which(beta!=0)
    Beta[[i]]=as.data.frame(as.matrix(beta))
    write.table(Beta[[i]],paste0('/home/unfated/simulation/coef_rid_trait',i,'.txt'))
  }
  snpid_rid=snpid
  Beta_rid = Beta
  write.table(as.character(snpid),'/home/unfated/simulation/snpid_rid.txt')
  write.table(as.character(Beta),'/home/unfated/simulation/Beta_rid.txt')
  
  #ridge + p value thresholding
  #gaussian
  Beta<-list()
  numselect<-NA
  k=20000
  snpid<-list()
  library(bigmemory)
  for(i in c(1:6,10)){
    i=5
    y2=YY[,i]
    # fit_cv = cv.glmnet(x1[1:k,],as.matrix(y2)[1:k,],family = 'gaussian',nfolds = 5,nlambda = 200,alpha = 0)
    # lambda0=fit_cv$lambda.1se
    # beta_cv=coef(fit_cv,s='lambda.1se')
    # print(sum(beta_cv!=0))
    # print(lambda0)
    lm<-lm(scale(y2)~1+scale(x1[,87]))
    summary(lm)
    lambda0=0.2
    fit_lar = glmnet(x1[1:10*k,],as.matrix(y2)[1:10*k,],lambda=lambda0,family = 'gaussian',alpha=0)
    beta = coef(fit_lar)
    # I=diag(x=lambda0,nrow=ncol(x1),ncol=ncol(x1))
    # sum(I)
    #x=big.matrix(x1,ncol=ncol(x1),nrow=nrow(x1))
    # R=t(x2)%*%x2
    # nrow(R)
    # V= chol2inv(chol(R))
    # H=x1%*%V%*%t(x1)
    # v=nrow(x1)-sum(diag(x1%*%solve(t(x1)%*%x1+I)%*%t(x1)))
    # chol2inv(chol(R))
    # v=nrow(x1)-sum(diag(2*H-H%*%t(H)))
    print(sum(beta!=0))
    print(sum(abs(beta)>0.01)-1)
    print(sum(abs(beta)>0.003333)-1)
    print(sum(abs(beta)>0.001)-1)
    id1=print(which(abs(beta)>0.003333))
    id2=print(which(abs(beta)>0.001))
    # id1=c(1, 65, 88, 255)
    # for(j in 2:length(id1)){
    # lm<-lm(scale(y2)~1+scale(x1[,256]))
    # print(summary(lm))
    # print(beta[id1[j]])
    # print(lm$coefficients)
    # }
    # for(i in 1:10){
    #   k= 0.000001*i
    # print(sum(beta^2>k)-1)
    # print(which(beta^2>k))
    # }

    #print(which(beta!=0))
    snpid[[i]]=which(beta!=0)
    Beta[[i]]=as.data.frame(as.matrix(beta))
    write.table(Beta[[i]],paste0('/home/unfated/simulation/coef_rid_trait',i,'.txt'))
  }
  snpid_rid=snpid
  Beta_rid = Beta
  write.table(as.character(snpid),'/home/unfated/simulation/snpid_rid.txt')
  write.table(as.character(Beta),'/home/unfated/simulation/Beta_rid.txt')
  
  # ridge p-value
  Beta<-list()
  numselect<-NA
  k=10000
  k2=c(1:1000)
  # k2=sample(1:10000,1000)
  # k2=c(1,10,100,1000)
  snpid<-list()
  library("lmridge") 
for(i in c(1:6,10)){
    i=1
    y2=YY[,i]
    fulldata=as.data.frame(x1)
    fulldata$y2<-y2
    ncol(fulldata)
    nrow(fulldata)
    fulldata1 = fulldata[,which(colSums(abs(fulldata)) !=0)]
    ncol(fulldata1)
    useddata=fulldata1[1:k,c(k2,ncol(fulldata1))]
    useddata = useddata[,which(colSums(abs(useddata)) !=0)]
    sum(is.na(useddata))
    ncol(useddata)
    nrow(useddata)
    ld=cor(useddata)-diag(1,ncol(useddata),ncol(useddata))
    length(which(ld>0.9))
    useddata1=useddata
    j=1
      while(j<ncol(useddata1)){
      ldd=cor(useddata1[,j],useddata1[,-j])
      if(length(which(ldd>0.9))>0){
      useddata1=useddata1[,-j]}else{j=j+1}
      }

    ncol(useddata1)
    nrow(useddata1)
useddata=useddata1[,-709]  
ncol(useddata)
nrow(useddata)

    # which(ld>0.99)
    useddata=scale(as.matrix(useddata))
    useddata=as.data.frame(useddata)
    sum(is.na(useddata))
    # which(is.na(useddata[2,]))
    # useddata$y2
  mod <- lmridge(y2 ~ ., data = useddata, scaling = "centered", K = 0.1) 
 
  Z=mod$Z
  Z=as.data.frame(Z)
  nrow(Z)
  Z=as.matrix(Z)
  nrow(Z)
  ncol(Z)
  X=as.matrix(useddata[,-ncol(X)])
  nrow(X)
  H=X%*%Z
  nrow(H)
  v=nrow(H)-sum(diag(H))
  # v
  ZZ=Z%*%t(Z)
  #ZZ
  Z%*%as.matrix(useddata[,ncol(useddata)])
  beta=as.vector(coef(mod))
  head(beta)
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
 qqplot(pval,runif(length(pval,0,1)))
 write.table(pval,'/home/unfated/pval.txt')
 num=length(which(pval<5*10^(-8)))
 num
min(pval) 
  # fit<-summary(mod) 
  # summary(mod)
  # str(fit)
  }