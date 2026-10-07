.libPaths(c("/home/user/rlib2", .libPaths())); suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
set.seed(3)
res <- NULL
for (r in 1:150) for (M in c(2e5)) {
  q <- 100
  blocks <- make_ld_ar1(rep(100,4), rho=0.6); m<-400; reg1<-1:200; reg2<-201:400
  w <- make_w(m, prop_nonzero=0.1, dist="lnorm"); Rb <- diag(m)
  gc <- make_gcov(q, type="cluster", rg=0.6, n_clusters=10); C <- make_intercept(q, overlap=1, rp=0.3*(gc$Rg>0))
  sim <- simulate_lrcpq(blocks, n=rep(3e5,q), w=w, Rb=Rb, gcov=gc$gcov, intercept=C, M=M)
  R1 <- as.matrix(bdiag(blocks[1:2])); R2 <- as.matrix(bdiag(blocks[3:4]))
  S1 <- which(w[reg1]>0); S2 <- which(w[reg2]>0)
  Z1 <- sim$Z[reg1,]; Z2 <- sim$Z[reg2,]
  f <- lrcp_distal(Z1,Z2,R1,R2,w[reg1],w[reg2],sim$n,gc$gcov,M,C,S1=S1,S2=S2,method="gls")
  g <- lrcp_gene(Z1[S1,],Z2[S2,],R1[S1,S1],R2[S2,S2],n=sim$n,gcov=gc$gcov,M=M,intercept=C,n_sim=0)
  g2 <- lrcp_gene(Z1[S1,],Z2[S2,],R1[S1,S1],R2[S2,S2],n=sim$n,gcov=gc$gcov,M=M,intercept=C,n_sim=0,
                  GA=(R1%*%(w[reg1]*R1))[S1,S1])
  zg<-as.vector(f$rho/f$se); zc<-as.vector(g$z); res <- rbind(res, data.frame(M=M, gls5=mean(abs(zg)>1.96,na.rm=T), cond5=mean(abs(zc)>1.96), gls3=mean(abs(zg)>qnorm(1-5e-4),na.rm=T), cond3=mean(abs(zc)>qnorm(1-5e-4)), gls4=mean(abs(zg)>qnorm(1-5e-5),na.rm=T), cond4=mean(abs(zc)>qnorm(1-5e-5)), n=length(zc)))
}
print(aggregate(.~M, res, mean)); saveRDS(res, commandArgs(TRUE)[1])
