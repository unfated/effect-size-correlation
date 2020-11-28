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
