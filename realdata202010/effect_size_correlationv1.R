#20201009#
library(tidyverse)
library(data.table)
GWAS.altas.genes <- fread("/home/unfated/realdata/202010gene/magma.P.r3.txt",stringsAsFactors = F)
#GWAS.altas.genes <- fread("C:/Users/22367/Desktop/fine mapping real data/magma.P.r3.txt",stringsAsFactors = F)
GWAS.altas.genes<-GWAS.altas.genes %>% as.data.frame()
a=ncol(GWAS.altas.genes)
colnames(GWAS.altas.genes)[2:a]<-paste0('GWAS',1:(a-1))
colnames(GWAS.altas.genes)[1:5]
colnames(GWAS.altas.genes)[1]<-'gene'
GWAS.altas.genes.all=GWAS.altas.genes
GWAS.altas.genes=GWAS.altas.genes[-1,-1]
nrow(GWAS.altas.genes)
GWAS.altas.genes.all[,1][1:100]
#GWAS.altas.genes=GWAS.altas.genes %>% t()
#col.names(GWAS.altas.genes)=GWAS.altas.genes.all$gene[]
# missing.gene=apply(GWAS.altas.genes,1,function(x)sum(is.na(x)))
# missing.gwas=apply(GWAS.altas.genes,2,function(x)sum(is.na(x)))
# head(missing.gene)

#----
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

## load the library 
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
## list the contents that are loaded into memory 
ls('package:TxDb.Hsapiens.UCSC.hg19.knownGene') 
## show the db object that is loaded by calling it's name 
TxDb.Hsapiens.UCSC.hg19.knownGene
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

#write.table(GWAS.altas.genes.all$gene[-1],'genename.txt',row.names=F,col.names=F,quote = F)

#map gene name----
library(org.Hs.eg.db)
ls('package:org.Hs.eg.db') 
org.Hs.egENSEMBL
## Bimap interface:
genename=GWAS.altas.genes.all$gene[-1]
x <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID 
mapped_genes <- mappedkeys(x) 
# Convert to a list 
ENS.name<- as.list(x[mapped_genes]) 
#if(length(xx) > 0) 
#{ # Get the Ensembl gene IDs for the first five genes 
 # xx[1:5] 
  # Get the first one 
  #xx[[1]] } 
ENS.name[1:5]

EG.name <- as.list(org.Hs.egENSEMBL2EG) 

EG.name[1:5]

EG.name['ENSG00000175899']
genename[1:5]
EG.name[genename[1:5]]
genename.EG<-EG.name[genename] 
head(genename.EG)
length(genename.EG)

genename.EG[[3]]

listtovec <- function(genename.EG) {
  genename.EG.ID<-NA
  for(i in 1:length(genename.EG)){
    a=genename.EG[[i]]
    if(length(a)==0){a=NA}
    genename.EG.ID[i]=a
  }
  return(genename.EG.ID)
}
genename.EG.ID=listtovec(genename.EG)
genename.EG.ID[1:100]

## Bimap interface: chromosome
x <- org.Hs.egCHR # Get the entrez gene identifiers that are mapped to a chromosome 
mapped_genes <- mappedkeys(x) 
# Convert to a list 
xx <- as.list(x[mapped_genes]) 
xx[1:5]
gene.chr=xx[genename.EG.ID]
gene.chr.ID=listtovec(gene.chr)
gene.chr.ID[20000:20100]
genename[1:5]
genename.EG[1:5]
genename.EG.ID[1:5]
gene.chr.ID[1:5]
gene.chr[1:5]
GWAS.altas.genes$chr=gene.chr.ID
rownames(GWAS.altas.genes)<-genename

#effect size correlation dist----
genes.chr1<-GWAS.altas.genes %>% subset(.,chr==1)
genes.chr1<-genes.chr1 %>% as.matrix() %>% t()
genes.chr1.1<-genes.chr1 

genes.chr1= data.frame(genes.chr1.1)
genes.chr2<-GWAS.altas.genes %>% subset(.,chr==2)
genes.chr2<-genes.chr2 %>% as.matrix() %>% t() %>%data.frame()

genes.chr1[1:3,1:3]
genes.chr2[1:3,1:3]
#genes.chr1 %>% as.data.frame.table() %>% 
#fwrite(.,'/home/unfated/realdata/202010gene/genes.chr1.txt',row.names=T,col.names=T)

write.table(genes.chr1,'/home/unfated/realdata/202010gene/genes.chr1.txt',row.names=T,col.names=T,quote=F)
#genes.chr2 %>% as.data.frame.table() %>% 
#fwrite(.,'/home/unfated/realdata/202010gene/genes.chr2.txt',row.names=T,col.names=T)
write.table(genes.chr2,'/home/unfated/realdata/202010gene/genes.chr2.txt',row.names=T,col.names=T,quote=F)


genes.chr1<-apply(genes.chr1,2,function(x) as.numeric(x))
genes.chr2<-apply(genes.chr2,2,function(x) as.numeric(x))

# new 20201112
library(tidyverse)
genes.chr1=read.table('/home/unfated/realdata/202010gene/genes.chr1.txt')
head(genes.chr1)
genes.chr2=read.table('/home/unfated/realdata/202010gene/genes.chr2.txt')
k=10
library(combinat, help, pos = 2, lib.loc = NULL)

count<-function(k) {
 # k=10
s1=sample(1:ncol(genes.chr1),k)
s1
s2=sample(1:ncol(genes.chr2),k)
s2
#-log((genes.chr1[,2]),10)%>%head()
#head(genes.chr1[,2])
#log(10,10)
genes.chr1.1<-apply(genes.chr1[,s1],2,function(x) -log(as.numeric(x),10))
genes.chr2.2<-apply(genes.chr2[,s2],2,function(x) -log(as.numeric(x),10))
genes.chr1.1[1:3,1:3]


corr.chr1.chr2=cor(genes.chr1.1,genes.chr2.2,use='pairwise.complete')
corr.chr1.chr2%>%nrow()
corr.chr1.chr2%>%ncol()
corr12=unlist(corr.chr1.chr2)
length(corr12)
head(corr12)
corr12[100]

#sum(corr12>0.01)
breakpoint=c(-1,-0.5,-0.25,-0.1,-0.05,-0.025,-0.01,-0.005,-0.001,0,0.001,0.005,0.01,0.025,0.05,0.1,0.25,0.5,1)
freq=cut(corr12,breaks=breakpoint)%>%table(.)
print('freq')
print(freq)
proportion=freq/length(corr12)
print('proportion')
print(proportion)
b=list(freq,proportion)
return(b)}

kk=c(10,20,50,100,200,300,400,500)

for(k in 1:10){
freq.n=list()
prop.n=list()
for(i in 1:length(kk)){
freq.n[[i]]=count(kk[i])[[1]]
prop.n[[i]]=count(kk[i])[[2]]
}


freq=freq.n%>%data.frame()
prop=prop.n%>%data.frame()
prop
freq.final=freq[,c(1,2*c(1:length(kk)))]
colnames(freq.final)=c('bin',kk)
prop.final=prop[,c(1,2*c(1:length(kk)))]
colnames(prop.final)=c('bin',kk)
write.csv(freq.final,paste0('/home/unfated/realdata/202010gene/freq.logcor20201113.real',k,'.csv'),row.names=F)
write.csv(prop.final,paste0('/home/unfated/realdata/202010gene/prop.logcor20201113.real',k,'.csv'),row.names=F)
print(k)
}



#sim----
n.size=ceiling(sum(!is.na(genes.chr1))/ncol(genes.chr1))

count2<-function(genes.chr1.1,genes.chr2.2){
  corr.chr1.chr2=cor(genes.chr1.1,genes.chr2.2,use='pairwise.complete')
corr12=unlist(corr.chr1.chr2)
 breakpoint=c(-1,-0.5,-0.25,-0.1,-0.05,-0.025,-0.01,-0.005,-0.001,0,0.001,0.005,0.01,0.025,0.05,0.1,0.25,0.5,1)
freq=cut(corr12,breaks=breakpoint)%>%table(.)
print('freq')
print(freq)
proportion=freq/length(corr12)
print('proportion')
print(proportion)
b=list(freq,proportion)
return(b)}

#simulation

 kk=c(10,20,50,100,200,300,400,500)

for(k in 1:10){
freq.sim=list()
prop.sim=list()
for(j in 1:length(kk)){
g1=matrix(0,nrow=n.size,ncol=kk[j])
g2=matrix(0,nrow=n.size,ncol=kk[j])
for(i in 1:kk[j]){
g1[,i]=-log(runif(n.size,-0,1),10)
g2[,i]=-log(runif(n.size,-0,1),10)}

freq.sim[[j]]=count2(g1,g2)[[1]]
prop.sim[[j]]=count2(g1,g2)[[2]]}

freq=freq.sim%>%data.frame()
prop=prop.sim%>%data.frame()
freq.final=freq[,c(1,2*c(1:length(kk)))]
colnames(freq.final)=c('bin',kk)
prop.final=prop[,c(1,2*c(1:length(kk)))]
colnames(prop.final)=c('bin',kk)
freq.final
prop.final
write.csv(freq.final,paste0('/home/unfated/realdata/202010gene/freq.logcor20201113.sim',k,'.csv'),row.names=F)
write.csv(prop.final,paste0('/home/unfated/realdata/202010gene/prop.logcor20201113.sim',k,'.csv'),row.names=F)
print(k)
}

#perm

count3<-function(k) {

 #k=20
s1=sample(1:ncol(genes.chr1),k)
s1
s2=sample(1:ncol(genes.chr2),k)
s2
#-log((genes.chr1[,2]),10)%>%head()
#head(genes.chr1[,2])
#log(10,10)
#k=10
a=apply(genes.chr1[,s1],2,function(x) -log(as.numeric(x),10))
b=apply(genes.chr2[,s2],2,function(x) -log(as.numeric(x),10))
resampling<-function(x) {
  aa=x
  aa[!is.na(x)]<-sample(na.omit(x),length(na.omit(x)),replace=F)
  return(aa)}
genes.chr1.1<-t(apply(a,1,resampling))
genes.chr2.2<-t(apply(b,1,resampling))



corr.chr1.chr2=cor(genes.chr1.1,genes.chr2.2,use='pairwise.complete')
corr.chr1.chr2%>%nrow()
corr.chr1.chr2%>%ncol()
corr12=unlist(corr.chr1.chr2)
length(corr12)

#sum(corr12>0.01)
breakpoint=c(-1,-0.5,-0.25,-0.1,-0.05,-0.025,-0.01,-0.005,-0.001,0,0.001,0.005,0.01,0.025,0.05,0.1,0.25,0.5,1)
freq=cut(corr12,breaks=breakpoint)%>%table(.)
print('freq')
print(freq)
proportion=freq/length(corr12)
print('proportion')
print(proportion)
b=list(freq,proportion)
return(b)}

#permu
kk=c(10,20,50,100,200,300,400,500)

for(k in 1:10){
freq.n=list()
prop.n=list()
for(i in 1:length(kk)){
freq.n[[i]]=count3(kk[i])[[1]]
prop.n[[i]]=count3(kk[i])[[2]]
}


freq=freq.n%>%data.frame()
prop=prop.n%>%data.frame()
prop
freq.final=freq[,c(1,2*c(1:length(kk)))]
colnames(freq.final)=c('bin',kk)
prop.final=prop[,c(1,2*c(1:length(kk)))]
colnames(prop.final)=c('bin',kk)
write.csv(freq.final,paste0('/home/unfated/realdata/202010gene/freq.logcor20201113.perm',k,'.csv'),row.names=F)
write.csv(prop.final,paste0('/home/unfated/realdata/202010gene/prop.logcor20201113.perm',k,'.csv'),row.names=F)
print(k)
}

