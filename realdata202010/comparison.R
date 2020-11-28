# # beta orginal
# list<-paste0('/home/unfated/Tim/beta.ls_fit',c(1:133),'.txt')
# BETA=c()
# for(i in 1:133){
# beta<-read.delim(list[i], stringsAsFactors=FALSE,sep = ' ')
# beta=beta[-1,]
# BETA=c(BETA,beta)}
# write.table(BETA,'/home/unfated/Tim/BETA.fit.txt')
# 
# #check NA snps
# IID2=list()
# IID3=c()
# trueb=read.delim('/home/unfated/Tim/effectsize_chr1.txt', stringsAsFactors=FALSE,sep = ' ')
# trueB=trueb[,1]
# for(j in 1:133){
#   data4<-NA
#   aa1=c('/home/unfated/chr1_block/block_')
#   bed1=paste0(aa1,j,'.bed')
#   bim1=paste0(aa1,j,'.bim')
#   fam1=paste0(aa1,j,'.fam')
#   block1=read.bed(bed1,bim1,fam1)
#   z_2=block1$snp[discovery,]
#   data4=data.frame(z_2)
#   id2=c()
#   for(l in 1:ncol(data4)){
#     if(sum(is.na(data4[,l]))>50000)
#       id2=c(id2,l)
#   }
#   IID2[[j]]=id2
# IID3[j]=ncol(data4)
# }
# 
# write.table(IID1,'remove_SNP_list.txt')
# 
# write.table(as.character(IID2),'remove_list_block.txt')
# length(IID1)
# 
# IID2=read.table('/home/unfated/remove_list_block.txt')
# IID2=read.table('C:/Users/22367/Desktop/fine mapping simulation/remove_list_block.txt',quote = "")
# dd=paste0('IID2[[',1:133,']]','=',as.character(IID2[1:133,]))
# write.table(dd,"IID2.txt",quote=F,row.names = F)
# 
# 
# # COMPARE
# 
# l2=c()
# for(i in 1:133){
# l2[i]=length(IID2[[i]])}
# 
# l1=c()
# for(i in 1:133){
#   beta<-read.delim(list[i], stringsAsFactors=FALSE,sep = ' ')
#   beta=beta[-1,]
#   l1[i]=length(beta)}
# l3=l1+l2
# IID3-l3
# sum(IID3)
# 
# # DUPLICATE
# list1<-paste0('C:/Users/22367/Desktop/fine mapping simulation/bim chr1/block_',c(1:133),'.bim')
# list1<-paste0('/home/unfated/chr1_block/block_',c(1:133),'.bim')
# Post=c()
# for(i in 1:133){
#   data<-read.delim(list1[i], stringsAsFactors=FALSE,quote = ' ',header = F)
#   dd=data$V4
#   Post=c(Post,dd)}

# data3<-read.table(file.choose(), stringsAsFactors=FALSE,header = T)
# nrow(data3)
# df=setdiff(Post,data3$POS)

# RECOVER na SNPS
# BETA=c()
# for(i in 1:133)
# {
#   beta<-read.delim(list[i], stringsAsFactors=FALSE,sep = ' ')
#   beta=beta[-1,]
#   if(length(IID2[[i]])>0)
#     for(j in 1:length(IID2[[i]])){
#       beta=append(beta,NA,after=IID2[[i]][j]-1)
#     }
#   BETA=c(BETA,beta)
# }
# #REMOVE DUPLICATE
# BETA=BETA[which(!duplicated(Post))]
# write.table(BETA,'/home/unfated/Tim/lasso_effectsize_chr1_59390SNP.txt')
# length(BETA)
# sum(is.na(BETA1))
# trueb=read.delim('/home/unfated/Tim/effectsize_chr1.txt', stringsAsFactors=FALSE,sep = ' ')
# trueB=trueb[,1]
# trueB=trueB[!is.na(BETA)]
# length(trueB)
# write.table(trueB,'/home/unfated/Tim/effectsize_chr1_remove_naSNP.txt')
# 
# trueb=read.table('/home/unfated/Tim/effectsize_chr1_remove_naSNP.txt')
# head(trueb)
# trueB<-trueb[,1]

holder=holder1
Rate <- function(trueB, holder,num) {
  list<-paste0(holder,c(1:133),'.txt')
  BETA=c()
  for(i in 1:num)
  {
    beta<-read.delim(list[i], stringsAsFactors=FALSE,sep = ' ')
    beta=beta[-1,]
    if(length(IID2[[i]])>0)
      for(j in 1:length(IID2[[i]])){
        beta=append(beta,NA,after=IID2[[i]][j]-1)
      }
    BETA=c(BETA,beta)
  }
  list1<-paste0('/home/unfated/chr1_block/block_',c(1:133),'.bim')
  Post=c()
  for(i in 1:num){
    data<-read.delim(list1[i], stringsAsFactors=FALSE,quote = ' ',header = F)
    dd=data$V4
    Post=c(Post,dd)}
  
  BETA=BETA[which(!duplicated(Post))]
  BETA=BETA[!is.na(BETA)]
  
  TN=0
  TP=0
  FP=0
  FN=0
  for(i in 1:length(BETA)){
    if(trueB[i]==0&& BETA[i]==0)
      TN=TN+1
    if(trueB[i]==0&& BETA[i]!=0)
      FP=FP+1
    if(trueB[i]!=0&& BETA[i]==0)
      FN=FN+1
    if(trueB[i]!=0&& BETA[i]!=0)
      TP=TP+1
  }
  print(paste0('FPR=',FP/(FP+TN)))
  print(paste0('FNR=',FN/(FN+TP)))
  print(paste0('FDR=',FP/(FP+TP)))

}

# blasso(X, y, T = 1000, thin = NULL, RJ = TRUE, M = NULL,
#        beta = NULL, lambda2 = 1, s2 = var(y-mean(y)),
#        case = c("default", "ridge", "hs", "ng"), mprior = 0, rd = NULL,
#        ab = NULL, theta = 0, rao.s2 = TRUE, icept = TRUE, 
#        normalize = TRUE, verb = 1)

