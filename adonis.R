
#-------------------------
#   健康与得病人群的adonis
#-------------------------
library(ggplot2)
library(vegan)
rm(list=ls())
result = rbind()
group = "group4"
sample = 'new_sample'
sample_map = read.table("../../00.data/mapping.file",sep="\t", header=T, check.names=F)
dt = read.table("../../00.data/new_mag.clust.profile.rename", sep="\t", header=T, check.names = F,row.names=1)

region = c("BJ", "SH")[2] # 在这边改变地域

sample_map = sample_map[which(sample_map$group2 == region),]
dtt = t(dt)
grps = unique(sample_map[,group])
com = t(combn(grps,2))

for (i in 1:nrow(com)){
  map = sample_map[sample_map[,group] %in% com[i,],]
  temp_dtt = dtt[map[,sample],]
  ado = adonis(temp_dtt~map[,group])
  r2 = ado$aov.tab$R2[1]
  p = ado$aov.tab$`Pr(>F)`[1]
  temp_result = data.frame(group = paste(com[i,1], com[i,2], sep=".vs." ),r2=r2,p=p, region = region)
  result = rbind(result, temp_result)
}

#write.table(result,"file.adonis.csv", sep=",",row.names=F)
#-----------上面运行完两次以后，再画图
result = result[order(result$region,result$r2, decreasing = T),]
result$group = factor(result$group, levels=rev(unique(result$group)))
result$sig = ifelse(result$p >= 0.05, "", ifelse( result$p < 0.05 & result$p >=0.01, "*", ifelse(result$p <0.01, "**", "***")))
ggplot(result,aes(y=group,x=r2, fill=region))+
  geom_bar(stat='identity',position=position_dodge2(width=1, preserve ='single' ))+
  geom_text(aes(label=sig),position=position_dodge2(width=1, preserve ='single' ), size=8, vjust=0.7, hjust=1)+
  theme_bw()+
  theme()

