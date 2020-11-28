set.seed(1)
r1=5
r2=4 #(h^2=0.75)
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
    # lambda=c(seq(0.01,lambdaused,length.out = 5)[-5],seq(lambdaused,2*lambdaused,length.out = 5))
    library(glmnet)
    fit<-cv.glmnet(as.matrix(useddata2),y,alpha = 0)
    # lambda=c(lambda,fit$lambda.min)
    lambda=c(seq(0.01,lambdaused,length.out = 5)[-5],seq(lambdaused,fit$lambda.min,by = lambdaused/4))
    lambda
    for(r3 in 1:length(lambda)){
    mod <- lmridge(y2 ~ ., data = useddata, scaling = "centered", K = lambda[r3]) 
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
    
    write.table(pval,paste0('/home/unfated/pval0816',herit,"_",r3-1,'.txt'))
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
    write.table(perform,paste0('/home/unfated/perform0816',herit,"_",r3-1,'.txt'))
    # perform[[k]]=c(FPR,FNR,FDR,power,herit/n.causal,n.causal,prop[k])
    # #length(which(BETA[,1]==0))
    # #length(which(BETA[,1]!=0))
    # print(paste0('k=',k))
    # 
    # store[[l]]=perform
    # # write.table(perform,paste0('/home/unfated/simulation/powerspectrum/SGFLasso-h',l,'.txt'))
    # 
  }