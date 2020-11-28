fig <- list()

library(ggplot2)
library(reshape2)
library(cowplot)
i=4
h2=c(0.1,0.25,0.5,0.75,0.9)
herit=h2[i]
lambdaused=(1-herit)/herit
lambdaused
lambda.min=11.77
lambda=c(seq(0.01,lambdaused,length.out = 5)[-5],seq(lambdaused,lambda.min,by = lambdaused/4))
lambda
# lambda=round(c(seq(0.01,lambdaused,length.out = 5)[-5],seq(lambdaused,2*lambdaused,length.out = 5)),17)
# lambda=c(lambda,11.77)
# for(j in 1:10){
#   prop = c(0.01,0.125,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95)
#   
#   m.causal=ceiling(883*prop)
#   perform <- read.csv(paste0("C:/Users/22367/Desktop/result/perform0816",h2[i],'_',j-1,".txt"), sep="")
#   perform$Fscore=2*perform$power*(1-perform$FPR)/(perform$power+(1-perform$FPR))
#   
#   perform$cutoff<-seq(1.5,6,by=0.5)
#   test_data_long1 <- reshape2::melt(perform[,c(1,2,3,4,5)], id="cutoff") 
#   p2<-ggplot(data=test_data_long1, aes(x=cutoff, y=value, colour=variable)) + geom_line()+geom_point()
#   p3<-p2+xlab('cutoff:-log10(pval)')+ggtitle(paste0('P=(h2=',h2[i],',lambda=',round(lambda[j],4), ',NO.causal=',m.causal[5],',perSNP h2=',round(h2[5]/m.causal[5],4),")"))
#   fig[[j]]=p3
# }
# 
# p5 <- cowplot::plot_grid(fig[[1]],fig[[2]],fig[[3]],fig[[4]],fig[[5]],fig[[6]],fig[[7]],fig[[8]],fig[[9]],fig[[10]],nrow = 3, labels = LETTERS[1:10])#将p1-p5四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D。（LETTERS[1:4] 意为提取26个大写英文字母的前四个:A、B、C、D）
# 
# p5
# png(width=1500,height=1000,paste0('ridge_cutoff0814lambda_',herit,'.png'))
# print(p5)
# dev.off()
# p5


#plot across lambda
fig<-list()
for(k in 1:10){
  perform1<-as.data.frame(matrix(1,length(lambda),5))
  
  
  for(j in 1:length(lambda)){
    prop = c(0.01,0.125,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95)
    
    m.causal=ceiling(883*prop)
    perform <- read.csv(paste0('C:/Users/22367/Desktop/result/perform0816',h2[i],'_',j-1,".txt"), sep="")
    perform$Fscore=2*perform$power*(1-perform$FPR)/(perform$power+(1-perform$FPR))
    
    perform$cutoff<-seq(1.5,6,by=0.5)
    perform1[j,]=perform[k,]
    # test_data_long1 <- reshape2::melt(perform[,c(1,2,3,4,5)], id="cutoff") 
    # p2<-ggplot(data=test_data_long1, aes(x=cutoff, y=value, colour=variable)) + geom_line()+geom_point()
    # p3<-p2+xlab('cutoff:-log10(pval)')+ggtitle(paste0('performance when P=(h2=',h2[5],',lambda=',round(lambda[j],4), ',NO.causal=',m.causal[5],',perSNP h2=',round(h2[5]/m.causal[5],4),")"))
    # fig[[j]]=p3
  }
  colnames(perform1)<-colnames(perform)
  perform1$lambda<-lambda
  test_data_long2 <- reshape2::melt(perform1[,c(1,2,3,4,6)], id="lambda") 
  p2<-ggplot(data=test_data_long2, aes(x=lambda, y=value, colour=variable)) + geom_line()+geom_point()+geom_vline(xintercept = lambdaused)+geom_vline(xintercept = lambda.min,color='red')
  p3<-p2+xlab('lambda')+ggtitle(paste0('P=[h2=',h2[i],',NO.causal=',m.causal[5],',cutoff(-log(p_val))=',perform$cutoff[k],"]"))
  fig[[k]]=p3
}

p5 <- cowplot::plot_grid(fig[[1]],fig[[2]],fig[[3]],fig[[4]],fig[[5]],fig[[6]],fig[[7]],fig[[8]],fig[[9]],fig[[10]], nrow = 5, labels = LETTERS[1:10])#将p1-p5四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D。（LETTERS[1:4] 意为提取26个大写英文字母的前四个:A、B、C、D）
# 
# p5
png(width=2000,height=1800,filename=paste0('ridge_0814_differlambdaCV_',herit,'.png'))
print(p5)
dev.off()
p5



