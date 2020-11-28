fig <- list()
i=10
for(j in 1:5){
prop = c(0.01,0.125,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95)
h2=c(0.1,0.25,0.5,0.75,0.9)
m.causal=ceiling(883*prop)
perform <- read.csv(paste0("C:/Users/22367/Desktop/perform",h2[j],'_',prop[i],".txt"), sep="")
perform$Fscore=2*perform$power*(1-perform$FPR)/(perform$power+(1-perform$FPR))
library(ggplot2)
library(reshape2)
perform$cutoff<-seq(1.5,6,by=0.5)
test_data_long1 <- reshape2::melt(perform[,c(1,2,3,4,5)], id="cutoff") 
p2<-ggplot(data=test_data_long1, aes(x=cutoff, y=value, colour=variable)) + geom_line()+geom_point()
p3<-p2+xlab('cutoff:-log10(pval)')+ggtitle(paste0('performance when P=(',h2[j],',',prop[i], ',',m.causal[i],',',round(h2[j]/m.causal[i],2),")"))
fig[[j]]=p3
}

p5 <- cowplot::plot_grid(fig[[1]],fig[[2]],fig[[3]],fig[[4]],fig[[5]], nrow = 3, labels = LETTERS[1:5])#将p1-p5四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D。（LETTERS[1:4] 意为提取26个大写英文字母的前四个:A、B、C、D）

png(paste0('ridge_cutoff',prop[i],'.png'))
print(p5)
dev.off()
p5
