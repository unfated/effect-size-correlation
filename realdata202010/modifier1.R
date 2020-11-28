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
bfile<-"/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train1"
#bfile <- useWTCCC(.bfile = F)

N <- nrow.bfile(bfile)
m <- 100000
discovery.sample <- sample(N, m)
# others.sample <- setdiff(1:N, discovery.sample)
# validation.sample <- sample(others.sample, validation.n)
# others.sample <- setdiff(1:N, c(discovery.sample, validation.sample))
# 
 discovery <- logical.vector(discovery.sample, N)
# validation <- logical.vector(validation.sample, N)

p <- ncol.bfile(bfile)
prop = 0.01
n.causal = ceiling(p*prop)
beta <- rep(0, p)
beta[sample(p, n.causal)] <- rnorm(n.causal)
write.table(beta,'effectsize_chr1.txt')
#Xb <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, keep=discovery)
# Xb <- pgs(bfile=bfile, weights=beta, chr=c(1,2), trace=1, keep=discovery)
Xb<- pgs(bfile=bfile, weights=beta, chr=c(1), trace=1, keep=discovery)
write.table(Xb,'pgs_chr1.txt')
herit <- 0.8
error <- rnorm(m, sd=sqrt(var(Xb)/herit*(1-herit)))
y <- Xb + error
write.table(y,'phenotype_chr1.txt')
Y <- rep(NA, N)

Y[discovery] <- y

#### Train ####
sum <- plink(bfile=bfile, cmd="--linear", plink2=TRUE, 
             pheno=Y, chr=c(1), silent=F, keep=discovery, 
            mem=fixed.mem)
summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)

est.beta <- summ$BETA
est.beta[is.na(est.beta)] <- 0
est.beta <- ifelse(summ$A1 == summ$ALT, est.beta, -est.beta)
write.table(est.beta,"est.beta.chr1.txt")
pval <- summ$P

pval[is.na(pval)] <- 1
write.table(pval,'pval.chr1.txt')

out = summ[,c(1,2,3,8,11)]
write.table(out,'summary.chr1.txt')



overlap <- (0:10)/10

pl <- pthresh.pipeline(beta=est.beta, pvals=pval, 
                       snp=summ$ID, keep.test = validation,
                       A1=summ$ALT, exclude.ambiguous = F, 
                       test.bfile = bfile, clump = TRUE, 
                       cluster=cl, trace=1, 
                       clump.options = list(r2=0.1))

Xb.validation <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, 
                     keep=validation)
error.validation <- rnorm(validation.n, sd=sqrt(var(Xb)/herit*(1-herit)))
y.validation <- Xb.validation + error.validation
v <- validate(pl, pheno=y.validation)
best.thres <- which(v$lambda == v$best.lambda)
best.beta <- v$best.beta

for(i in 1:length(overlap)) {
  # i <- 1
  test.sample <- integer(0)
  overlap.n <- test.n * overlap[i]
  others.n <- test.n  - overlap.n
  test.sample <- c(sample(discovery.sample, overlap.n), 
                   sample(others.sample, others.n))
  test <- logical.vector(test.sample, N)
  
  zero <- pl$beta[[1]]$beta == 0
  pval2 <- pval
  pval2[zero] <- 1
  
  o <- rank(pval2)
  est.beta100 <- est.beta
  est.beta100[o > 100] <- 0
  BETA <- cbind(beta, best.beta, est.beta100, est.beta)
  
  PGS.test <- pgs(bfile=bfile, weights=BETA, cluster=cl, keep=test)
  
  cor.test <- as.vector(cor(PGS.test)[1,-1])
  
  fill.in.results.table(n=test.n, seed=attr(filename, "seed"),
                        m=m, n.causal=n.causal, 
                        cor = cor.test, 
                        type2 = c("best", "100", "full"), 
                        p = colSums(BETA[,-1] != 0), 
                        P = nrow(BETA[,-1]), 
                        overlap=overlap[i])

}


save(results.table, file=paste0(filename, ".RData"))
