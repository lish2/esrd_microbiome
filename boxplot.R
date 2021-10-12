dt = read.table("../../../00.data/new_mag.clust.profile.rename",sep="\t", header=T,check.names = F, row.names = 1)
map = read.table("../../../00.data/mapping.file",sep="\t", header=T, check.names=F)
sm = read.table("../../../00.data/strain.taxo.txt", sep="\t", header=T, check.names = F)

map = map[which(map$group2 != "BJ"), c(5,2)]
dt = dt[,map$new_sample]

sm = sm[which(sm$family %in% c("Lachnospiraceae", "Oscillospiraceae", "Ruminococcaceae", "Prevotellaceae")),]
sm = sm[,c(1,5)]
library(reshape2)

dt$strain = rownames(dt)
dl = melt(dt)
data = merge(dl, sm, by.x='strain', by.y='user_genome')
data2 = merge(data, map, by.x='variable', by.y='sample')
library(ggplot2)

data2$value = log10(data2$value+1)
ggplot(data2, aes(x=group1, y=value, fill=family))+
  geom_boxplot()+
  theme_bw()+
  ylab("log10 abundance")




#-------------------------------------------
#             上海的
#-------------------------------------------
rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggsignif)

dt = read.table("../../../00.data/new_mag.clust.profile.rename",sep="\t", header=T,check.names = F, row.names = 1)
map = read.table("../../../00.data/mapping.file",sep="\t", header=T, check.names=F)
sm = read.table("../../../00.data/strain.taxo.txt", sep="\t", header=T, check.names = F)

map = map[which(map$group2 != "BJ"), c(7,5)]
dt = dt[,map$new_sample]


sm = sm[which(sm[,'enriched'] !="na" ),]
sm = sm[,c(1,9)]

#dt = scale(dt)# scaled
dt = dt[sm$user_genome,]
dt = data.frame(dt)
dsm = aggregate(. ~ sm$enriched, dt, sum)

dl = melt(dsm)
data = merge(dl, map, by.x='variable', by.y='new_sample')

comp = combn(unique(data[, "group4"]),2,list)
data$group4 = factor(data$group4, levels=c("Control", "CKD","CKD-5nD","ESRD"))
#data$value = log10(data$value)

sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

data1 = data[which(data$`sm$enriched` == "ESRD-enriched"),]
data2 = data[which(data$`sm$enriched` == "HC-enriched"),]


p1 = ggplot(data1, aes(x=group4, y=value))+
  geom_boxplot(fill='red')+
  theme_bw()+
  geom_signif(comparisons = comp,test='wilcox.test',step_increase = 0.1)


p2 = ggplot(data2, aes(x=group4, y=value))+
  geom_boxplot(fill='blue')+
  theme_bw()+
  geom_signif(comparisons = comp,test='wilcox.test',step_increase = 0.1)

library(patchwork)

p1+p2

