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
p
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
library('lsgl')

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
sum(is.na(YY))
for(j in 1:traitnum){
  md=median(YY2[,j],na.rm=T)
  YY2[is.na(YY2[,j]),j]=md}
YY3=YY2+2
# marginal GWAS
for(j in 5:10){
  Y=YY2[,j]
  print(mean(Y))}
numGWAS<-NA
#case rate
for(j in 5:10){
  Y=YY3[,j]
  sum <- plink(bfile=bfile, cmd="--glm", plink2=T, 
               pheno=Y, silent=F, keep=discovery
              )
  summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)
  SNP<-summ[,c('ID','P')]
  pval <- summ$P
  pval[is.na(pval)] <- 1
  #clumping
  clump(pval,summ$POS,r2=0.5)
  aaa=0
  for(i in 1:length(pval)){
    if(pval[i]<5*10^(-8)){aaa=aaa+1}
  }
  aaa
  numGWAS[j]=aaa
  print(aaa)
}
numGWAS

#pure Multivariate Lasso
y1=as.matrix(YY[1:5])
lambda<-lsgl::lambda(x1,y1, alpha=0, lambda.min=0.05)
fit<-lsgl::fit(x = x1, y = y1, alpha = 0, lambda = lambda)
betalasso=as.matrix(coef(fit, traitnum)[,-1])
