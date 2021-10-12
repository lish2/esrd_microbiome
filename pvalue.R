
#-------------------------
#  分组的菌株含量差异
#-------------------------
library(ggplot2)
library(vegan)
rm(list=ls())
result = rbind()
sample_map = read.table("../../00.data/mapping.file",sep="\t", header=T, check.names=F)
dt = read.table("../00.data/mag.clust.profile", sep="\t", header=T, check.names = F,row.names=1)


region = c("BJ", "SH")[2] # 在这边改变地域

sample_map = sample_map[which(sample_map$protopathy !="na" & sample_map$group2 == region & sample_map$group3 != "CKD"),]


dt = dt[,sample_map$sample]
grps = unique(sample_map$protopathy)
com = t(combn(grps,2))
nspecies = nrow(dt)
names = rownames(dt)
for (n in 1:nspecies){
  temp_dt = dt[n,]
  for(c in 1:nrow(com)){
    g1 = com[c,1]
    g2 = com[c,2]
    g1s = sample_map[which(sample_map$protopathy == g1),1]
    g2s = sample_map[which(sample_map$protopathy == g2),1]
    dt1 = as.matrix(temp_dt[,g1s])
    dt2 = as.matrix(temp_dt[,g2s])
    m1 = mean(dt1)
    m2 = mean(dt2)
    p = wilcox.test(dt1,dt2)$p.value
    temp_result = data.frame(name = names[n],g1=g1, g2=g2,mean1 = m1, mean2=m2,pvalue=p)
    result = rbind(result, temp_result)
  }
}

write.table(result, paste(region, "pvalue.csv", sep=""), sep=",", row.names=F)

#write.table(fdrtool(as.matrix(read.table("clipboard", sep="\t"))[,1], stat='pvalue')$qval, "clipboard-128",sep="\t", row.names=F,col.names = F)




#-------------------------
#  上海 分组的菌株含量差异
#-------------------------
library(ggplot2)
library(vegan)
rm(list=ls())

sample_map = read.table("../../00.data/mapping.file",sep="\t", header=T, check.names=F)
dt = read.table("../../00.data/new_mag.clust.profile.rename", sep="\t", header=T, check.names = F,row.names=1)
dt = dt/colSums(dt)*100

sample_map = sample_map[which(sample_map$group2 %in% c("SH")),]

# 哪个分组进行比较
group = 'group4'
dt = dt[,sample_map$new_sample]
grps = unique(sample_map[,group])
com = t(combn(grps,2))
nspecies = nrow(dt)
names = rownames(dt)
result = rbind()
for (n in 1:nspecies){
  temp_dt = dt[n,]
  for(c in 1:nrow(com)){
    g1 = com[c,1]
    g2 = com[c,2]
    g1s = sample_map[which(sample_map[,group] == g1),'new_sample']
    g2s = sample_map[which(sample_map[,group] == g2),'new_sample']
    dt1 = as.matrix(temp_dt[,g1s])
    dt2 = as.matrix(temp_dt[,g2s])
    c1 = sum(dt1 != 0 )
    c2 = sum(dt2 != 0)
    m1 = mean(dt1)
    m2 = mean(dt2)
    p = wilcox.test(dt1,dt2)$p.value
    temp_result = data.frame(name = names[n],g1=g1, g2=g2,mean1 = m1, mean2=m2,pvalue=p, count1=c1, count2= c2)
    result = rbind(result, temp_result)
  }
}


write.table(result, paste("SH", ".new_sample_group4.pvalue.csv", sep=""), sep=",", row.names=F)
