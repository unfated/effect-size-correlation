#  qsub -I -q psychipc -l nodes=1:ppn=4,walltime=15:00:00,mem=20gb  
a0='#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=1gb
#PBS -l walltime=60:00:00
#PBS -m abe
#PBS -q small_ext
#PBS -N block
'
a='module load plink
cd /psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/
'
b='
plink --bfile  whiteBritish1_train1 --chr 1 --from-bp '
d=' --to-bp ' 
f=' --make-bed --out /home/unfated/chr1_block/block_'
data <- read.delim("C:/Users/22367/Desktop/fine mapping simulation/fourier_ls-chr1.bed", stringsAsFactors=FALSE)
data2=data[,-1]
tt=c()
for(i in 1:nrow(data2)){
tt=paste0(tt,b,data2[i,1],d,data2[i,2],f,i)
}
write.table(paste0(a0,a,tt),'block.pbs',row.names=F,col.names =F,quote=F)
