

itr=c(1:1)
option=c('lasso')
holder=c('/home/unfated/Tim/')
input=c('/home/unfated/chr1_block/block_')
impute=2
y=read.table('/home/unfated/Tim/phenotype_chr1.txt')
y=y$V1
length(y)

Lasso <- function(holder,itr,option,discovery,impute,y) {
  library('KRIS')
  library('glmnet')
  library('monomvn')
  IID=list()
  for(i in itr){
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
    if(sum(is.na(data1[,k]))>50000)
      id=c(id,k)
  }
  IID[[i]]=id
  data3=data1
  if(length(id>0))
    data3=data1[,-id]
  for(j in 1:ncol(data3)){
    data3[is.na(data3[,j]),j]=1
  }
  x1=as.matrix(data3)
  cv1 = cv.glmnet(x1, y)
  if(option=='lasso'||option=='both'){
          beta.ls_fit1=coef(cv1, s = "lambda.1se")
         write.table(as.matrix(beta.ls_fit1),paste0(holder,'beta.ls_fit',i,'.txt'))}
  if(option=='blasso'||option=='both' ){
          bls<-blasso(x1, y, T = 500, thin = NULL, RJ = TRUE, M = NULL,
         beta = NULL, lambda2 = 1, s2 = var(y-mean(y)),
         case = c("default"), mprior = 0, rd = NULL,
         ab = NULL, theta = 0, rao.s2 = TRUE, icept = TRUE, 
         normalize = TRUE, verb = 1)
          
  write.table(as.matrix(beta.ls_fit1),paste0(holder,'beta.bls_fit',i,'.txt'))}
  }}
  